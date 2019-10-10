% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
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
	% StartTime: 2019-10-10 09:55:40
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:40
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(9);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0, 0; -t44, -t45, 0, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0, 0; t45, t44, 0, 0, 0, 0; t43, t46, 0, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:41
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (101->28), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t256 = cos(qJ(4));
	t257 = cos(qJ(1));
	t276 = t256 * t257;
	t255 = sin(qJ(1));
	t275 = qJD(1) * t255;
	t274 = qJD(1) * t257;
	t254 = sin(qJ(4));
	t273 = qJD(2) * t254;
	t272 = qJD(2) * t255;
	t271 = qJD(2) * t256;
	t270 = qJD(2) * t257;
	t269 = qJD(4) * t254;
	t268 = qJD(4) * t256;
	t267 = qJD(4) * t257;
	t253 = qJ(2) + pkin(9);
	t251 = sin(t253);
	t266 = t251 * t271;
	t265 = t251 * t272;
	t252 = cos(t253);
	t264 = t252 * t272;
	t263 = t251 * t270;
	t262 = t252 * t270;
	t261 = qJD(4) * t252 - qJD(1);
	t260 = qJD(1) * t252 - qJD(4);
	t259 = t261 * t254;
	t258 = t260 * t255 + t263;
	t250 = -t260 * t276 + (t259 + t266) * t255;
	t249 = t261 * t256 * t255 + (t260 * t257 - t265) * t254;
	t248 = t258 * t256 + t257 * t259;
	t247 = t258 * t254 - t261 * t276;
	t1 = [t250, -t256 * t262 + (t254 * t267 + t256 * t275) * t251, 0, t247, 0, 0; -t248, -t256 * t264 + (t255 * t269 - t256 * t274) * t251, 0, -t249, 0, 0; 0, -t252 * t269 - t266, 0, -t251 * t268 - t252 * t273, 0, 0; t249, t254 * t262 + (-t254 * t275 + t256 * t267) * t251, 0, t248, 0, 0; t247, t254 * t264 + (t254 * t274 + t255 * t268) * t251, 0, t250, 0, 0; 0, t251 * t273 - t252 * t268, 0, t251 * t269 - t252 * t271, 0, 0; -t251 * t274 - t264, -t252 * t275 - t263, 0, 0, 0, 0; -t251 * t275 + t262, t252 * t274 - t265, 0, 0, 0, 0; 0, qJD(2) * t252, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:41
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (161->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->35)
	t280 = cos(qJ(1));
	t278 = qJ(2) + pkin(9);
	t276 = cos(t278);
	t292 = qJD(4) * t276;
	t286 = -qJD(1) + t292;
	t302 = t280 * t286;
	t285 = qJD(1) * t276 - qJD(4);
	t274 = sin(t278);
	t279 = sin(qJ(1));
	t296 = qJD(2) * t279;
	t290 = t274 * t296;
	t301 = t285 * t280 - t290;
	t300 = qJD(1) * t279;
	t299 = qJD(1) * t280;
	t298 = qJD(2) * t274;
	t297 = qJD(2) * t276;
	t295 = qJD(2) * t280;
	t277 = qJ(4) + pkin(10);
	t273 = sin(t277);
	t294 = qJD(4) * t273;
	t293 = qJD(4) * t274;
	t275 = cos(t277);
	t291 = t275 * t293;
	t289 = t276 * t296;
	t288 = t274 * t295;
	t287 = t276 * t295;
	t284 = t286 * t279;
	t283 = t274 * t299 + t289;
	t282 = -t274 * t300 + t287;
	t281 = t285 * t279 + t288;
	t272 = t273 * t284 - t301 * t275;
	t271 = t301 * t273 + t275 * t284;
	t270 = t273 * t302 + t281 * t275;
	t269 = t281 * t273 - t275 * t302;
	t1 = [t272, -t275 * t287 + (t275 * t300 + t280 * t294) * t274, 0, t269, 0, 0; -t270, -t275 * t289 + (-t275 * t299 + t279 * t294) * t274, 0, -t271, 0, 0; 0, -t273 * t292 - t275 * t298, 0, -t273 * t297 - t291, 0, 0; t271, t282 * t273 + t280 * t291, 0, t270, 0, 0; t269, t283 * t273 + t279 * t291, 0, t272, 0, 0; 0, t273 * t298 - t275 * t292, 0, t273 * t293 - t275 * t297, 0, 0; -t283, -t276 * t300 - t288, 0, 0, 0, 0; t282, t276 * t299 - t290, 0, 0, 0, 0; 0, t297, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:42
	% EndTime: 2019-10-10 09:55:42
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (161->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->35)
	t341 = sin(qJ(1));
	t340 = qJ(2) + pkin(9);
	t338 = cos(t340);
	t347 = qJD(1) * t338 - qJD(4);
	t336 = sin(t340);
	t342 = cos(qJ(1));
	t357 = qJD(2) * t342;
	t350 = t336 * t357;
	t364 = t347 * t341 + t350;
	t358 = qJD(2) * t341;
	t352 = t336 * t358;
	t363 = t347 * t342 - t352;
	t362 = qJD(1) * t341;
	t361 = qJD(1) * t342;
	t360 = qJD(2) * t336;
	t359 = qJD(2) * t338;
	t339 = qJ(4) + pkin(10);
	t335 = sin(t339);
	t356 = qJD(4) * t335;
	t355 = qJD(4) * t336;
	t354 = qJD(4) * t338;
	t337 = cos(t339);
	t353 = t337 * t355;
	t351 = t338 * t358;
	t349 = t338 * t357;
	t348 = qJD(1) - t354;
	t346 = t348 * t341;
	t345 = t348 * t342;
	t344 = -t336 * t361 - t351;
	t343 = t336 * t362 - t349;
	t334 = t335 * t346 + t363 * t337;
	t333 = -t363 * t335 + t337 * t346;
	t332 = t335 * t345 - t364 * t337;
	t331 = t364 * t335 + t337 * t345;
	t1 = [-t334, -t337 * t349 + (t337 * t362 + t342 * t356) * t336, 0, t331, 0, 0; t332, -t337 * t351 + (-t337 * t361 + t341 * t356) * t336, 0, t333, 0, 0; 0, -t335 * t354 - t337 * t360, 0, -t335 * t359 - t353, 0, 0; t344, -t338 * t362 - t350, 0, 0, 0, 0; -t343, t338 * t361 - t352, 0, 0, 0, 0; 0, t359, 0, 0, 0, 0; t333, t343 * t335 - t342 * t353, 0, t332, 0, 0; -t331, t344 * t335 - t341 * t353, 0, t334, 0, 0; 0, -t335 * t360 + t337 * t354, 0, -t335 * t355 + t337 * t359, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end