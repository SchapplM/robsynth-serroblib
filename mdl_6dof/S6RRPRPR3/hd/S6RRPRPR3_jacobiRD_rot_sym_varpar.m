% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:07:54
	% EndTime: 2019-10-10 10:07:54
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
	t253 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:07:54
	% EndTime: 2019-10-10 10:07:55
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (161->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->35)
	t280 = cos(qJ(1));
	t278 = qJ(2) + pkin(10);
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
	t277 = qJ(4) + pkin(11);
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
	% StartTime: 2019-10-10 10:07:54
	% EndTime: 2019-10-10 10:07:55
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (337->31), mult. (233->49), div. (0->0), fcn. (233->6), ass. (0->36)
	t326 = qJ(2) + pkin(10);
	t322 = sin(t326);
	t325 = qJD(4) + qJD(6);
	t348 = t322 * t325;
	t323 = cos(t326);
	t347 = t323 * t325;
	t327 = sin(qJ(1));
	t346 = qJD(1) * t327;
	t328 = cos(qJ(1));
	t345 = qJD(1) * t328;
	t344 = qJD(2) * t322;
	t343 = qJD(2) * t323;
	t342 = qJD(2) * t327;
	t341 = qJD(2) * t328;
	t324 = qJ(4) + pkin(11) + qJ(6);
	t320 = sin(t324);
	t340 = t320 * t348;
	t321 = cos(t324);
	t339 = t321 * t348;
	t338 = t327 * t325 * t321;
	t337 = t328 * t325 * t320;
	t336 = t322 * t342;
	t335 = t322 * t341;
	t334 = -qJD(1) + t347;
	t333 = qJD(1) * t323 - t325;
	t332 = t320 * t334;
	t331 = t322 * t345 + t323 * t342;
	t330 = t322 * t346 - t323 * t341;
	t329 = t333 * t327 + t335;
	t316 = -t321 * t343 + t340;
	t315 = -t320 * t343 - t339;
	t314 = t327 * t332 + (-t333 * t328 + t336) * t321;
	t313 = -t320 * t336 - t337 - t321 * t346 + (t320 * t345 + t338) * t323;
	t312 = t329 * t321 + t328 * t332;
	t311 = -t334 * t328 * t321 + t329 * t320;
	t1 = [t314, t330 * t321 + t322 * t337, 0, t311, 0, t311; -t312, -t331 * t321 + t327 * t340, 0, -t313, 0, -t313; 0, -t320 * t347 - t321 * t344, 0, t315, 0, t315; t313, -t330 * t320 + t328 * t339, 0, t312, 0, t312; t311, t331 * t320 + t322 * t338, 0, t314, 0, t314; 0, t320 * t344 - t321 * t347, 0, t316, 0, t316; -t331, -t323 * t346 - t335, 0, 0, 0, 0; -t330, t323 * t345 - t336, 0, 0, 0, 0; 0, t343, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end