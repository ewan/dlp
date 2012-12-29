package org.jqgibbs;



public class GenericSampler extends Sampler {
	
	public GenericSampler(Model m) {
		super(m);
	}
	
	@Override
	public ChainLink variate() {
		ChainLink newLink = new ChainLink();
		for (RandomVar <?> gv : this.current) {
			newLink.add(gv.cloneShallow());
		}
		for (RandomVar <?> gv : newLink) {
			gv.updatePosterior(newLink);
		}
		this.setChainLink(newLink);
		return newLink;
	}
}
