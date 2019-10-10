% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t35 = cos(pkin(9));
	t34 = sin(pkin(9));
	t33 = (-t34 * t37 + t35 * t36) * qJD(1);
	t32 = (t34 * t36 + t35 * t37) * qJD(1);
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t32, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0; -t32, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:56
	% EndTime: 2019-10-09 23:30:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (27->9), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t71 = sin(pkin(9));
	t72 = cos(pkin(9));
	t74 = sin(qJ(1));
	t76 = cos(qJ(1));
	t67 = t76 * t71 - t74 * t72;
	t73 = sin(qJ(5));
	t80 = qJD(5) * t73;
	t75 = cos(qJ(5));
	t79 = qJD(5) * t75;
	t66 = -t74 * t71 - t76 * t72;
	t65 = t67 * qJD(1);
	t78 = -t65 * t73 + t66 * t79;
	t77 = t65 * t75 + t66 * t80;
	t64 = t66 * qJD(1);
	t63 = t64 * t75 - t67 * t80;
	t62 = t64 * t73 + t67 * t79;
	t1 = [t78, 0, 0, 0, t63, 0; t62, 0, 0, 0, t77, 0; 0, 0, 0, 0, t79, 0; -t77, 0, 0, 0, -t62, 0; t63, 0, 0, 0, t78, 0; 0, 0, 0, 0, -t80, 0; t64, 0, 0, 0, 0, 0; t65, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:57
	% EndTime: 2019-10-09 23:30:57
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (108->26), mult. (317->48), div. (0->0), fcn. (353->8), ass. (0->31)
	t298 = sin(pkin(9));
	t299 = cos(pkin(9));
	t300 = sin(qJ(1));
	t301 = cos(qJ(1));
	t272 = t301 * t298 - t300 * t299;
	t277 = sin(qJ(6));
	t279 = cos(qJ(6));
	t270 = t272 * qJD(1);
	t271 = -t300 * t298 - t301 * t299;
	t278 = sin(qJ(5));
	t280 = cos(qJ(5));
	t294 = qJD(5) * t280;
	t286 = t270 * t278 - t271 * t294;
	t282 = qJD(6) * t272 - t286;
	t269 = t271 * qJD(1);
	t293 = qJD(6) * t278;
	t289 = t271 * t293 - t269;
	t302 = t282 * t277 + t289 * t279;
	t297 = t277 * t280;
	t296 = t279 * t280;
	t295 = qJD(5) * t278;
	t292 = qJD(6) * t280;
	t288 = -t272 * t293 + t270;
	t287 = t269 * t278 + t272 * t294;
	t285 = t277 * t292 + t279 * t295;
	t284 = t277 * t295 - t279 * t292;
	t283 = -qJD(6) * t271 + t287;
	t281 = -t289 * t277 + t282 * t279;
	t268 = t288 * t277 + t283 * t279;
	t267 = -t283 * t277 + t288 * t279;
	t1 = [t281, 0, 0, 0, t269 * t296 - t285 * t272, t267; t268, 0, 0, 0, t270 * t296 + t285 * t271, t302; 0, 0, 0, 0, -t277 * t293 + t279 * t294, -t284; -t302, 0, 0, 0, -t269 * t297 + t284 * t272, -t268; t267, 0, 0, 0, -t270 * t297 - t284 * t271, t281; 0, 0, 0, 0, -t277 * t294 - t279 * t293, -t285; t270 * t280 + t271 * t295, 0, 0, 0, t287, 0; -t269 * t280 + t272 * t295, 0, 0, 0, t286, 0; 0, 0, 0, 0, t295, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end