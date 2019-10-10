% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR11_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(6));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(6));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t122 = sin(qJ(2));
	t123 = sin(qJ(1));
	t130 = t122 * t123;
	t125 = cos(qJ(1));
	t129 = t122 * t125;
	t124 = cos(qJ(2));
	t128 = t123 * t124;
	t127 = t124 * t125;
	t120 = sin(pkin(6));
	t126 = qJD(1) * t120;
	t121 = cos(pkin(6));
	t1 = [0, t125 * t126, (-t121 * t130 + t127) * qJD(2) + (t121 * t127 - t130) * qJD(1), 0, 0, 0; 0, t123 * t126, (t121 * t129 + t128) * qJD(2) + (t121 * t128 + t129) * qJD(1), 0, 0, 0; 0, 0, t120 * qJD(2) * t122, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->9), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t147 = sin(qJ(2));
	t148 = sin(qJ(1));
	t156 = t147 * t148;
	t150 = cos(qJ(1));
	t155 = t147 * t150;
	t149 = cos(qJ(2));
	t154 = t148 * t149;
	t153 = t149 * t150;
	t145 = sin(pkin(6));
	t152 = qJD(1) * t145;
	t151 = t145 * qJD(2) * t147;
	t146 = cos(pkin(6));
	t144 = (t146 * t155 + t154) * qJD(2) + (t146 * t154 + t155) * qJD(1);
	t143 = (t146 * t156 - t153) * qJD(2) + (-t146 * t153 + t156) * qJD(1);
	t1 = [0, t150 * t152, -t143, 0, t143, 0; 0, t148 * t152, t144, 0, -t144, 0; 0, 0, t151, 0, -t151, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:20
	% EndTime: 2019-10-10 12:10:20
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (64->31), mult. (200->61), div. (0->0), fcn. (214->10), ass. (0->32)
	t298 = qJD(5) - qJD(3);
	t264 = sin(pkin(6));
	t267 = sin(qJ(3));
	t294 = t264 * t267;
	t271 = cos(qJ(3));
	t293 = t264 * t271;
	t273 = cos(qJ(1));
	t292 = t264 * t273;
	t268 = sin(qJ(2));
	t269 = sin(qJ(1));
	t291 = t268 * t269;
	t290 = t268 * t273;
	t272 = cos(qJ(2));
	t289 = t269 * t272;
	t288 = t273 * t272;
	t287 = qJD(1) * t264;
	t286 = qJD(2) * t264;
	t285 = t269 * t287;
	t284 = t273 * t287;
	t283 = t268 * t286;
	t265 = cos(pkin(6));
	t279 = t265 * t288 - t291;
	t278 = t265 * t289 + t290;
	t262 = t265 * t290 + t289;
	t277 = t265 * t291 - t288;
	t270 = cos(qJ(5));
	t266 = sin(qJ(5));
	t261 = -t277 * qJD(1) + t279 * qJD(2);
	t260 = t278 * qJD(1) + t262 * qJD(2);
	t259 = -t262 * qJD(1) - t278 * qJD(2);
	t258 = -t279 * qJD(1) + t277 * qJD(2);
	t1 = [0, t284, -t258, 0, t258, (t259 * t271 + t267 * t284) * t266 - (t259 * t267 - t271 * t284) * t270 - t298 * ((t267 * t277 + t269 * t293) * t266 - (t269 * t294 - t271 * t277) * t270); 0, t285, t260, 0, -t260, (t261 * t271 + t267 * t285) * t266 - (t261 * t267 - t271 * t285) * t270 + t298 * ((t262 * t267 + t271 * t292) * t266 + (t262 * t271 - t267 * t292) * t270); 0, 0, t283, 0, -t283, (t271 * t266 - t267 * t270) * t272 * t286 + t298 * ((-t265 * t271 + t268 * t294) * t266 + (t265 * t267 + t268 * t293) * t270);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end