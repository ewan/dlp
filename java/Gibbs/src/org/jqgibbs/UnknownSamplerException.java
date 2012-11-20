package org.jqgibbs;

public class UnknownSamplerException extends GibbsException {

	private static final long serialVersionUID = 8246639443259717460L;

	public UnknownSamplerException() {
		super();
	}

	public UnknownSamplerException(String message) {
		super(message);
	}

	public UnknownSamplerException(Throwable cause) {
		super(cause);
	}

	public UnknownSamplerException(String message, Throwable cause) {
		super(message, cause);
	}

}
