% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:02:28
	% EndTime: 2019-10-10 09:02:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:02:28
	% EndTime: 2019-10-10 09:02:28
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
	% StartTime: 2019-10-10 09:02:28
	% EndTime: 2019-10-10 09:02:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(11);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:02:28
	% EndTime: 2019-10-10 09:02:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(11);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0, 0; -t35, 0, -t36, 0, 0, 0; 0, 0, -t44, 0, 0, 0; t36, 0, t35, 0, 0, 0; t34, 0, t37, 0, 0, 0; 0, 0, -t43, 0, 0, 0; -qJD(1) * t38, 0, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:02:30
	% EndTime: 2019-10-10 09:02:30
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (108->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t245 = cos(qJ(4));
	t246 = cos(qJ(3));
	t261 = qJD(1) * t246;
	t251 = -qJD(4) + t261;
	t264 = t245 * t251;
	t257 = qJD(4) * t246;
	t252 = -qJD(1) + t257;
	t243 = sin(qJ(4));
	t244 = sin(qJ(3));
	t260 = qJD(3) * t244;
	t254 = t243 * t260;
	t263 = t252 * t245 - t254;
	t262 = qJD(1) * t244;
	t259 = qJD(3) * t246;
	t258 = qJD(4) * t244;
	t256 = t243 * t262;
	t255 = t245 * t262;
	t253 = t245 * t260;
	t250 = t251 * t243;
	t249 = t243 * t258 - t245 * t259;
	t248 = t243 * t259 + t245 * t258;
	t247 = t252 * t243 + t253;
	t242 = qJ(1) + pkin(11);
	t241 = cos(t242);
	t240 = sin(t242);
	t239 = t247 * t240 - t241 * t264;
	t238 = t263 * t240 + t241 * t250;
	t237 = t240 * t264 + t247 * t241;
	t236 = t240 * t250 - t263 * t241;
	t1 = [t239, 0, t240 * t255 + t249 * t241, t236, 0, 0; -t237, 0, t249 * t240 - t241 * t255, -t238, 0, 0; 0, 0, -t243 * t257 - t253, -t248, 0, 0; t238, 0, -t240 * t256 + t248 * t241, t237, 0, 0; t236, 0, t248 * t240 + t241 * t256, t239, 0, 0; 0, 0, -t245 * t257 + t254, t249, 0, 0; -t240 * t259 - t241 * t262, 0, -t240 * t261 - t241 * t260, 0, 0, 0; -t240 * t262 + t241 * t259, 0, -t240 * t260 + t241 * t261, 0, 0, 0; 0, 0, t259, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:02:30
	% EndTime: 2019-10-10 09:02:30
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (264->31), mult. (233->50), div. (0->0), fcn. (233->6), ass. (0->33)
	t300 = qJ(1) + pkin(11);
	t295 = sin(t300);
	t301 = qJ(4) + qJ(5);
	t298 = cos(t301);
	t321 = t295 * t298;
	t296 = cos(t300);
	t297 = sin(t301);
	t320 = t296 * t297;
	t299 = qJD(4) + qJD(5);
	t302 = sin(qJ(3));
	t319 = t299 * t302;
	t303 = cos(qJ(3));
	t318 = t299 * t303;
	t317 = qJD(1) * t302;
	t316 = qJD(1) * t303;
	t315 = qJD(3) * t302;
	t314 = qJD(3) * t303;
	t313 = t297 * t317;
	t312 = t298 * t317;
	t311 = t297 * t315;
	t310 = t298 * t315;
	t309 = t295 * t315;
	t308 = -qJD(1) + t318;
	t307 = -t299 + t316;
	t306 = t295 * t307;
	t291 = t297 * t319 - t298 * t314;
	t305 = t297 * t314 + t298 * t319;
	t304 = t308 * t297 + t310;
	t289 = -t307 * t298 * t296 + t304 * t295;
	t288 = -t297 * t309 - t299 * t320 - qJD(1) * t321 + (qJD(1) * t320 + t299 * t321) * t303;
	t287 = t304 * t296 + t298 * t306;
	t286 = t297 * t306 + (-t308 * t298 + t311) * t296;
	t1 = [t289, 0, t291 * t296 + t295 * t312, t286, t286, 0; -t287, 0, t291 * t295 - t296 * t312, -t288, -t288, 0; 0, 0, -t297 * t318 - t310, -t305, -t305, 0; t288, 0, -t295 * t313 + t305 * t296, t287, t287, 0; t286, 0, t305 * t295 + t296 * t313, t289, t289, 0; 0, 0, -t298 * t318 + t311, t291, t291, 0; -t295 * t314 - t296 * t317, 0, -t295 * t316 - t296 * t315, 0, 0, 0; -t295 * t317 + t296 * t314, 0, t296 * t316 - t309, 0, 0, 0; 0, 0, t314, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:02:30
	% EndTime: 2019-10-10 09:02:30
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (488->32), mult. (293->48), div. (0->0), fcn. (293->6), ass. (0->34)
	t316 = qJ(4) + qJ(5) + qJ(6);
	t311 = sin(t316);
	t317 = qJ(1) + pkin(11);
	t314 = cos(t317);
	t337 = t311 * t314;
	t312 = cos(t316);
	t313 = sin(t317);
	t336 = t312 * t313;
	t315 = qJD(4) + qJD(5) + qJD(6);
	t318 = sin(qJ(3));
	t335 = t315 * t318;
	t319 = cos(qJ(3));
	t334 = t315 * t319;
	t333 = qJD(1) * t318;
	t332 = qJD(1) * t319;
	t331 = qJD(3) * t318;
	t330 = qJD(3) * t319;
	t329 = t311 * t335;
	t328 = t312 * t335;
	t327 = t313 * t331;
	t326 = t314 * t331;
	t325 = -qJD(1) + t334;
	t324 = -t315 + t332;
	t323 = t311 * t325;
	t322 = t313 * t330 + t314 * t333;
	t321 = t313 * t333 - t314 * t330;
	t320 = t324 * t313 + t326;
	t307 = -t312 * t330 + t329;
	t306 = -t311 * t330 - t328;
	t305 = t313 * t323 + (-t324 * t314 + t327) * t312;
	t304 = -t315 * t337 - t311 * t327 - qJD(1) * t336 + (qJD(1) * t337 + t315 * t336) * t319;
	t303 = t320 * t312 + t314 * t323;
	t302 = -t325 * t314 * t312 + t320 * t311;
	t1 = [t305, 0, t321 * t312 + t314 * t329, t302, t302, t302; -t303, 0, -t322 * t312 + t313 * t329, -t304, -t304, -t304; 0, 0, -t311 * t334 - t312 * t331, t306, t306, t306; t304, 0, -t321 * t311 + t314 * t328, t303, t303, t303; t302, 0, t322 * t311 + t313 * t328, t305, t305, t305; 0, 0, t311 * t331 - t312 * t334, t307, t307, t307; -t322, 0, -t313 * t332 - t326, 0, 0, 0; -t321, 0, t314 * t332 - t327, 0, 0, 0; 0, 0, t330, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end