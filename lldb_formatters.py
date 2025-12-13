"""
LLDB Pretty Printer for math++ Matrix class

This module provides custom formatters for the Matrix template class
to display matrices in a human-readable format during debugging sessions.

Usage:
    In LLDB, run: command script import /path/to/lldb_formatters.py
    Or add to .lldbinit file for automatic loading
"""

try:
    import lldb
except ImportError:
    # LLDB module not available when running standalone
    pass

import re


def _unwrap_swig(obj):
    """Unwrap SWIG proxy objects to get the actual LLDB object."""
    # Try multiple unwrapping strategies
    if hasattr(obj, '__dict__') and 'this' in obj.__dict__:
        return obj.__dict__['this']
    if hasattr(obj, 'this'):
        return getattr(obj, 'this')
    # Already unwrapped or different type
    return obj


def matrix_summary(valobj, internal_dict):
    """
    Provides a summary string for Matrix objects.
    
    Args:
        valobj: LLDB SBValue object representing the Matrix
        internal_dict: Internal dictionary for caching
        
    Returns:
        A summary string showing matrix dimensions and type
    """
    try:
        # Try to get the actual LLDB SBValue - handle SWIG proxy
        val = valobj
        
        # Check if it's callable (sometimes SWIG objects are)
        if callable(valobj):
            val = valobj()
        
        # Try direct unwrap
        if not hasattr(val, 'GetTypeName'):
            val = _unwrap_swig(valobj)
        
        # Last resort: try to import lldb and cast
        if not hasattr(val, 'GetTypeName'):
            try:
                import lldb
                # The valobj might need to be casted
                if str(type(valobj)) == "<class 'SwigPyObject'>":
                    return "Matrix<?, ?> (SWIG binding issue - use 'p' instead of 'fr v')"
            except:
                pass
        
        # Get type name
        type_name = val.GetTypeName()
        
        # Parse Matrix<COLUMNS, ROWS, T> pattern
        match = re.match(r'Matrix<(\d+),\s*(\d+),\s*(.+)>', type_name)
        if not match:
            return f"Matrix ({type_name})"
        
        columns = int(match.group(1))
        rows = int(match.group(2))
        value_type = match.group(3).strip()
        
        # Get the data array
        data_member = val.GetChildMemberWithName('data')
        
        if not data_member or not data_member.IsValid():
            return f"Matrix<{columns}x{rows}, {value_type}>"
        
        # For small matrices, show a preview of values
        if columns <= 4 and rows <= 4:
            try:
                values = []
                for col in range(columns):
                    col_array = data_member.GetChildAtIndex(col)
                    if col_array and col_array.IsValid():
                        for row in range(rows):
                            val_elem = col_array.GetChildAtIndex(row)
                            if val_elem and val_elem.IsValid():
                                val_str = val_elem.GetValue()
                                if val_str:
                                    values.append(val_str)
                
                if len(values) == columns * rows:
                    preview = ', '.join(values[:min(6, len(values))])
                    suffix = '...' if len(values) > 6 else ''
                    return f"Matrix<{columns}x{rows}, {value_type}> [{preview}{suffix}]"
            except:
                pass
        
        return f"Matrix<{columns}x{rows}, {value_type}>"
        
    except Exception as e:
        return f"Matrix<?, ?> (error: {str(e)})"


class MatrixSyntheticProvider:
    """
    Synthetic children provider for Matrix template class.
    
    This class allows LLDB to display matrix elements in a structured way,
    showing individual elements with their row/column positions.
    """
    
    def __init__(self, valobj, internal_dict):
        """
        Initialize the synthetic provider.
        
        Args:
            valobj: LLDB SBValue object representing the Matrix
            internal_dict: Internal dictionary for caching
        """
        # Unwrap SWIG proxy
        self.valobj = _unwrap_swig(valobj)
        self.update()
    
    def num_children(self):
        """
        Returns the number of children (matrix elements).
        
        Returns:
            Total number of elements (columns * rows)
        """
        return self.columns * self.rows
    
    def get_child_index(self, name):
        """
        Get the index for a child by name.
        
        Args:
            name: Child name (e.g., "[0,0]", "[1,2]")
            
        Returns:
            Index of the child or -1 if not found
        """
        try:
            # Parse names like "[col,row]" or "data"
            if name == "data":
                return -1
            
            match = re.match(r'\[(\d+),(\d+)\]', name)
            if match:
                col = int(match.group(1))
                row = int(match.group(2))
                if 0 <= col < self.columns and 0 <= row < self.rows:
                    return row * self.columns + col
        except:
            pass
        return -1
    
    def get_child_at_index(self, index):
        """
        Get the child element at the given index.
        
        Args:
            index: Linear index (0 to columns*rows-1)
            
        Returns:
            LLDB SBValue representing the matrix element
        """
        if index < 0 or index >= self.num_children():
            return None
        
        try:
            # Convert linear index to column/row
            col = index % self.columns
            row = index // self.columns
            
            # Access data[col][row]
            col_array = self.data_member.GetChildAtIndex(col)
            if not col_array.IsValid():
                return None
            
            element = col_array.GetChildAtIndex(row)
            if not element.IsValid():
                return None
            
            # Create a named value for display
            name = f"[{col},{row}]"
            return element.CreateChildAtOffset(name, 0, element.GetType())
            
        except Exception as e:
            return None
    
    def update(self):
        """
        Update internal state when the variable changes.
        
        This method is called by LLDB to refresh the provider's state.
        """
        try:
            # Extract template parameters from type name
            type_name = self.valobj.GetTypeName()
            match = re.match(r'Matrix<(\d+),\s*(\d+),\s*(.+)>', type_name)
            
            if match:
                self.columns = int(match.group(1))
                self.rows = int(match.group(2))
                self.value_type = match.group(3).strip()
            else:
                self.columns = 0
                self.rows = 0
                self.value_type = "unknown"
            
            # Get the data member
            self.data_member = self.valobj.GetChildMemberWithName('data')
            
        except Exception as e:
            self.columns = 0
            self.rows = 0
            self.value_type = "error"
            self.data_member = None
    
    def has_children(self):
        """
        Indicates whether this object has children.
        
        Returns:
            True if the matrix has elements
        """
        return self.num_children() > 0


# If running standalone, provide usage information
if __name__ == "__main__":
    print(__doc__)
    print("\nThis script should be imported into LLDB, not run directly.")
    print("Usage: (lldb) command script import /path/to/lldb_formatters.py")
