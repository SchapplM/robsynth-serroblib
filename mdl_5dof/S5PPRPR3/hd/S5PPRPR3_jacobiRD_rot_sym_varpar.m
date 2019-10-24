% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->15), div. (0->0), fcn. (24->6), ass. (0->11)
	t66 = sin(pkin(7));
	t69 = sin(qJ(3));
	t75 = t66 * t69;
	t70 = cos(qJ(3));
	t74 = t66 * t70;
	t68 = cos(pkin(7));
	t73 = t68 * t69;
	t72 = t68 * t70;
	t71 = qJD(3) * sin(pkin(8));
	t67 = cos(pkin(8));
	t1 = [0, 0, (-t67 * t72 - t75) * qJD(3), 0, 0; 0, 0, (-t67 * t74 + t73) * qJD(3), 0, 0; 0, 0, -t70 * t71, 0, 0; 0, 0, (t67 * t73 - t74) * qJD(3), 0, 0; 0, 0, (t67 * t75 + t72) * qJD(3), 0, 0; 0, 0, t69 * t71, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->10)
	t76 = sin(pkin(7));
	t77 = cos(pkin(8));
	t81 = t76 * t77;
	t78 = cos(pkin(7));
	t80 = t77 * t78;
	t79 = qJD(3) * sin(pkin(8));
	t74 = qJ(3) + pkin(9);
	t73 = cos(t74);
	t72 = sin(t74);
	t1 = [0, 0, (-t72 * t76 - t73 * t80) * qJD(3), 0, 0; 0, 0, (t72 * t78 - t73 * t81) * qJD(3), 0, 0; 0, 0, -t73 * t79, 0, 0; 0, 0, (t72 * t80 - t73 * t76) * qJD(3), 0, 0; 0, 0, (t72 * t81 + t73 * t78) * qJD(3), 0, 0; 0, 0, t72 * t79, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:51
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->24), mult. (140->57), div. (0->0), fcn. (148->8), ass. (0->27)
	t227 = sin(pkin(8));
	t231 = sin(qJ(5));
	t241 = t227 * t231;
	t232 = cos(qJ(5));
	t240 = t227 * t232;
	t228 = sin(pkin(7));
	t229 = cos(pkin(8));
	t239 = t228 * t229;
	t226 = qJ(3) + pkin(9);
	t224 = sin(t226);
	t230 = cos(pkin(7));
	t238 = t230 * t224;
	t225 = cos(t226);
	t237 = t230 * t225;
	t236 = qJD(3) * t225;
	t235 = qJD(5) * t231;
	t234 = qJD(5) * t232;
	t233 = t227 * qJD(3) * t224;
	t223 = t228 * t224 + t229 * t237;
	t221 = t225 * t239 - t238;
	t222 = t228 * t225 - t229 * t238;
	t220 = -t224 * t239 - t237;
	t219 = t223 * qJD(3);
	t218 = t222 * qJD(3);
	t217 = t221 * qJD(3);
	t216 = t220 * qJD(3);
	t1 = [0, 0, -t219 * t232 - t222 * t235, 0, -t218 * t231 + (-t223 * t232 - t230 * t241) * qJD(5); 0, 0, -t217 * t232 - t220 * t235, 0, -t216 * t231 + (-t221 * t232 - t228 * t241) * qJD(5); 0, 0, (t224 * t235 - t232 * t236) * t227, 0, t231 * t233 + (-t225 * t240 + t229 * t231) * qJD(5); 0, 0, t219 * t231 - t222 * t234, 0, -t218 * t232 + (t223 * t231 - t230 * t240) * qJD(5); 0, 0, t217 * t231 - t220 * t234, 0, -t216 * t232 + (t221 * t231 - t228 * t240) * qJD(5); 0, 0, (t224 * t234 + t231 * t236) * t227, 0, t232 * t233 + (t225 * t241 + t229 * t232) * qJD(5); 0, 0, t218, 0, 0; 0, 0, t216, 0, 0; 0, 0, -t233, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end