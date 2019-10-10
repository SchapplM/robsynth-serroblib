% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t14 = qJD(1) * sin(qJ(1));
	t13 = qJD(1) * cos(qJ(1));
	t1 = [-t13, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = cos(pkin(9));
	t61 = sin(pkin(9));
	t60 = (t61 * t64 - t62 * t63) * qJD(1);
	t59 = (t61 * t63 + t62 * t64) * qJD(1);
	t1 = [-t59, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->10), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t67 = sin(qJ(4));
	t74 = qJD(4) * t67;
	t69 = cos(qJ(4));
	t73 = qJD(4) * t69;
	t65 = sin(pkin(9));
	t66 = cos(pkin(9));
	t68 = sin(qJ(1));
	t70 = cos(qJ(1));
	t62 = t70 * t65 - t68 * t66;
	t63 = t68 * t65 + t70 * t66;
	t60 = t63 * qJD(1);
	t72 = -t60 * t67 + t62 * t73;
	t71 = -t60 * t69 - t62 * t74;
	t61 = t62 * qJD(1);
	t59 = t61 * t69 - t63 * t74;
	t58 = -t61 * t67 - t63 * t73;
	t1 = [t71, 0, 0, t58, 0, 0; t59, 0, 0, t72, 0, 0; 0, 0, 0, t74, 0, 0; -t72, 0, 0, -t59, 0, 0; t58, 0, 0, t71, 0, 0; 0, 0, 0, t73, 0, 0; -t61, 0, 0, 0, 0, 0; -t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:24
	% EndTime: 2019-10-09 23:39:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (44->11), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->18)
	t82 = qJ(4) + pkin(10);
	t80 = sin(t82);
	t90 = qJD(4) * t80;
	t81 = cos(t82);
	t89 = qJD(4) * t81;
	t83 = sin(pkin(9));
	t84 = cos(pkin(9));
	t85 = sin(qJ(1));
	t86 = cos(qJ(1));
	t77 = t86 * t83 - t85 * t84;
	t78 = t85 * t83 + t86 * t84;
	t75 = t78 * qJD(1);
	t88 = -t75 * t80 + t77 * t89;
	t87 = -t75 * t81 - t77 * t90;
	t76 = t77 * qJD(1);
	t74 = t76 * t81 - t78 * t90;
	t73 = -t76 * t80 - t78 * t89;
	t1 = [t87, 0, 0, t73, 0, 0; t74, 0, 0, t88, 0, 0; 0, 0, 0, t90, 0, 0; -t88, 0, 0, -t74, 0, 0; t73, 0, 0, t87, 0, 0; 0, 0, 0, t89, 0, 0; -t76, 0, 0, 0, 0, 0; -t75, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:39:26
	% EndTime: 2019-10-09 23:39:26
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (162->26), mult. (317->51), div. (0->0), fcn. (353->8), ass. (0->35)
	t321 = sin(pkin(9));
	t322 = cos(pkin(9));
	t323 = sin(qJ(1));
	t324 = cos(qJ(1));
	t291 = t324 * t321 - t323 * t322;
	t299 = sin(qJ(6));
	t300 = cos(qJ(6));
	t290 = -t323 * t321 - t324 * t322;
	t288 = t290 * qJD(1);
	t298 = qJ(4) + pkin(10);
	t297 = cos(t298);
	t296 = sin(t298);
	t318 = qJD(4) * t296;
	t307 = -t288 * t297 + t291 * t318;
	t303 = -qJD(6) * t290 + t307;
	t289 = t291 * qJD(1);
	t314 = qJD(6) * t297;
	t308 = t291 * t314 + t289;
	t325 = t303 * t299 - t308 * t300;
	t320 = t296 * t299;
	t319 = t296 * t300;
	t317 = qJD(4) * t297;
	t316 = qJD(4) * t299;
	t315 = qJD(4) * t300;
	t313 = qJD(6) * t299;
	t312 = qJD(6) * t300;
	t309 = t290 * t314 + t288;
	t306 = t289 * t297 + t290 * t318;
	t305 = -t296 * t313 + t297 * t315;
	t304 = t296 * t312 + t297 * t316;
	t302 = qJD(6) * t291 + t306;
	t301 = -t308 * t299 - t303 * t300;
	t287 = t309 * t299 + t302 * t300;
	t286 = -t302 * t299 + t309 * t300;
	t1 = [t301, 0, 0, -t289 * t319 + t305 * t290, 0, t286; t287, 0, 0, t288 * t319 + t305 * t291, 0, -t325; 0, 0, 0, t296 * t315 + t297 * t313, 0, t304; t325, 0, 0, t289 * t320 - t304 * t290, 0, -t287; t286, 0, 0, -t288 * t320 - t304 * t291, 0, t301; 0, 0, 0, -t296 * t316 + t297 * t312, 0, t305; t288 * t296 + t291 * t317, 0, 0, t306, 0, 0; t289 * t296 - t290 * t317, 0, 0, t307, 0, 0; 0, 0, 0, -t317, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end