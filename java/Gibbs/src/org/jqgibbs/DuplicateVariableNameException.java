package org.jqgibbs;

public class DuplicateVariableNameException extends GibbsException {

	private static final long serialVersionUID = 4847602360055140030L;

	public DuplicateVariableNameException() {
		super();
	}

	public DuplicateVariableNameException(String message) {
		super(message);
	}

	public DuplicateVariableNameException(Throwable cause) {
		super(cause);
	}

	public DuplicateVariableNameException(String message, Throwable cause) {
		super(message, cause);
	}

}
