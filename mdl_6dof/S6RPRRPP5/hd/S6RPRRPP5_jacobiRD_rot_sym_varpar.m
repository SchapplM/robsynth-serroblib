% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPP5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
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
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(9));
	t17 = sin(pkin(9));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:40
	% EndTime: 2019-10-10 01:16:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(9) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0, 0; -t37, 0, -t38, 0, 0, 0; 0, 0, -qJD(3) * t40, 0, 0, 0; t38, 0, t37, 0, 0, 0; t36, 0, t39, 0, 0, 0; 0, 0, -qJD(3) * t41, 0, 0, 0; -t48, 0, 0, 0, 0, 0; t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:42
	% EndTime: 2019-10-10 01:16:42
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (101->28), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t245 = cos(qJ(4));
	t246 = cos(qJ(1));
	t265 = t245 * t246;
	t244 = sin(qJ(1));
	t264 = qJD(1) * t244;
	t263 = qJD(1) * t246;
	t243 = sin(qJ(4));
	t262 = qJD(3) * t243;
	t261 = qJD(3) * t244;
	t260 = qJD(3) * t245;
	t259 = qJD(3) * t246;
	t258 = qJD(4) * t243;
	t257 = qJD(4) * t245;
	t256 = qJD(4) * t246;
	t242 = pkin(9) + qJ(3);
	t240 = sin(t242);
	t255 = t240 * t260;
	t254 = t240 * t261;
	t241 = cos(t242);
	t253 = t241 * t261;
	t252 = t240 * t259;
	t251 = t241 * t259;
	t250 = qJD(4) * t241 - qJD(1);
	t249 = qJD(1) * t241 - qJD(4);
	t248 = t250 * t243;
	t247 = t249 * t244 + t252;
	t239 = -t249 * t265 + (t248 + t255) * t244;
	t238 = t250 * t245 * t244 + (t249 * t246 - t254) * t243;
	t237 = t247 * t245 + t246 * t248;
	t236 = t247 * t243 - t250 * t265;
	t1 = [t239, 0, -t245 * t251 + (t243 * t256 + t245 * t264) * t240, t236, 0, 0; -t237, 0, -t245 * t253 + (t244 * t258 - t245 * t263) * t240, -t238, 0, 0; 0, 0, -t241 * t258 - t255, -t240 * t257 - t241 * t262, 0, 0; t238, 0, t243 * t251 + (-t243 * t264 + t245 * t256) * t240, t237, 0, 0; t236, 0, t243 * t253 + (t243 * t263 + t244 * t257) * t240, t239, 0, 0; 0, 0, t240 * t262 - t241 * t257, t240 * t258 - t241 * t260, 0, 0; -t240 * t263 - t253, 0, -t241 * t264 - t252, 0, 0, 0; -t240 * t264 + t251, 0, t241 * t263 - t254, 0, 0, 0; 0, 0, qJD(3) * t241, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:42
	% EndTime: 2019-10-10 01:16:42
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (101->28), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t309 = cos(qJ(1));
	t305 = pkin(9) + qJ(3);
	t304 = cos(t305);
	t311 = qJD(1) * t304 - qJD(4);
	t328 = t311 * t309;
	t307 = sin(qJ(1));
	t303 = sin(t305);
	t321 = qJD(3) * t309;
	t314 = t303 * t321;
	t327 = t311 * t307 + t314;
	t326 = qJD(1) * t307;
	t325 = qJD(1) * t309;
	t306 = sin(qJ(4));
	t324 = qJD(3) * t306;
	t323 = qJD(3) * t307;
	t308 = cos(qJ(4));
	t322 = qJD(3) * t308;
	t320 = qJD(4) * t306;
	t319 = qJD(4) * t308;
	t318 = qJD(4) * t309;
	t317 = t303 * t322;
	t316 = t303 * t323;
	t315 = t304 * t323;
	t313 = t304 * t321;
	t312 = -qJD(4) * t304 + qJD(1);
	t310 = t312 * t309;
	t302 = t308 * t328 + (t312 * t306 - t317) * t307;
	t301 = t312 * t308 * t307 + (t316 - t328) * t306;
	t300 = t306 * t310 - t327 * t308;
	t299 = t327 * t306 + t308 * t310;
	t1 = [-t302, 0, -t308 * t313 + (t306 * t318 + t308 * t326) * t303, t299, 0, 0; t300, 0, -t308 * t315 + (t307 * t320 - t308 * t325) * t303, t301, 0, 0; 0, 0, -t304 * t320 - t317, -t303 * t319 - t304 * t324, 0, 0; -t303 * t325 - t315, 0, -t304 * t326 - t314, 0, 0, 0; -t303 * t326 + t313, 0, t304 * t325 - t316, 0, 0, 0; 0, 0, qJD(3) * t304, 0, 0, 0; t301, 0, -t306 * t313 + (t306 * t326 - t308 * t318) * t303, t300, 0, 0; -t299, 0, -t306 * t315 + (-t306 * t325 - t307 * t319) * t303, t302, 0, 0; 0, 0, -t303 * t324 + t304 * t319, -t303 * t320 + t304 * t322, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:16:42
	% EndTime: 2019-10-10 01:16:42
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t276 = cos(qJ(1));
	t272 = pkin(9) + qJ(3);
	t271 = cos(t272);
	t278 = qJD(1) * t271 - qJD(4);
	t295 = t278 * t276;
	t274 = sin(qJ(1));
	t270 = sin(t272);
	t288 = qJD(3) * t276;
	t281 = t270 * t288;
	t294 = t278 * t274 + t281;
	t293 = qJD(1) * t274;
	t292 = qJD(1) * t276;
	t273 = sin(qJ(4));
	t291 = qJD(3) * t273;
	t290 = qJD(3) * t274;
	t275 = cos(qJ(4));
	t289 = qJD(3) * t275;
	t287 = qJD(4) * t273;
	t286 = qJD(4) * t275;
	t285 = qJD(4) * t276;
	t284 = t270 * t289;
	t283 = t270 * t290;
	t282 = t271 * t290;
	t280 = t271 * t288;
	t279 = -qJD(4) * t271 + qJD(1);
	t277 = t279 * t276;
	t269 = t275 * t295 + (t279 * t273 - t284) * t274;
	t268 = t279 * t275 * t274 + (t283 - t295) * t273;
	t267 = t273 * t277 - t294 * t275;
	t266 = t294 * t273 + t275 * t277;
	t1 = [-t269, 0, -t275 * t280 + (t273 * t285 + t275 * t293) * t270, t266, 0, 0; t267, 0, -t275 * t282 + (t274 * t287 - t275 * t292) * t270, t268, 0, 0; 0, 0, -t271 * t287 - t284, -t270 * t286 - t271 * t291, 0, 0; t268, 0, -t273 * t280 + (t273 * t293 - t275 * t285) * t270, t267, 0, 0; -t266, 0, -t273 * t282 + (-t273 * t292 - t274 * t286) * t270, t269, 0, 0; 0, 0, -t270 * t291 + t271 * t286, -t270 * t287 + t271 * t289, 0, 0; t270 * t292 + t282, 0, t271 * t293 + t281, 0, 0, 0; t270 * t293 - t280, 0, -t271 * t292 + t283, 0, 0, 0; 0, 0, -qJD(3) * t271, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end