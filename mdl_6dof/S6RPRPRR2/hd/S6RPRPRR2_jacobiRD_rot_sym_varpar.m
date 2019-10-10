% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
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
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(10);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(10);
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
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t57 = qJ(1) + pkin(10);
	t53 = sin(t57);
	t62 = qJD(1) * t53;
	t55 = cos(t57);
	t61 = qJD(1) * t55;
	t56 = qJ(3) + pkin(11);
	t52 = sin(t56);
	t60 = qJD(3) * t52;
	t54 = cos(t56);
	t59 = qJD(3) * t54;
	t58 = qJD(3) * t55;
	t51 = t53 * t60 - t54 * t61;
	t50 = t52 * t61 + t53 * t59;
	t49 = t52 * t58 + t54 * t62;
	t48 = t52 * t62 - t54 * t58;
	t1 = [t51, 0, t48, 0, 0, 0; -t49, 0, -t50, 0, 0, 0; 0, 0, -t60, 0, 0, 0; t50, 0, t49, 0, 0, 0; t48, 0, t51, 0, 0, 0; 0, 0, -t59, 0, 0, 0; -t62, 0, 0, 0, 0, 0; t61, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:48:00
	% EndTime: 2019-10-10 00:48:00
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (161->28), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->34)
	t270 = cos(qJ(5));
	t267 = qJ(3) + pkin(11);
	t265 = cos(t267);
	t288 = qJD(1) * t265;
	t273 = -qJD(5) + t288;
	t291 = t270 * t273;
	t274 = qJD(5) * t265 - qJD(1);
	t263 = sin(t267);
	t269 = sin(qJ(5));
	t282 = qJD(3) * t269;
	t278 = t263 * t282;
	t290 = t274 * t270 - t278;
	t268 = qJ(1) + pkin(10);
	t264 = sin(t268);
	t289 = t263 * t264;
	t266 = cos(t268);
	t287 = qJD(1) * t266;
	t286 = qJD(1) * t269;
	t285 = qJD(1) * t270;
	t284 = qJD(3) * t265;
	t283 = qJD(3) * t266;
	t281 = qJD(3) * t270;
	t280 = qJD(5) * t269;
	t279 = qJD(5) * t270;
	t277 = t263 * t281;
	t276 = t264 * t284;
	t275 = t265 * t283;
	t272 = t273 * t269;
	t271 = t274 * t269 + t277;
	t262 = t271 * t264 - t266 * t291;
	t261 = t290 * t264 + t266 * t272;
	t260 = t264 * t291 + t271 * t266;
	t259 = t264 * t272 - t290 * t266;
	t1 = [t262, 0, -t270 * t275 + (t264 * t285 + t266 * t280) * t263, 0, t259, 0; -t260, 0, -t270 * t276 + (t264 * t280 - t266 * t285) * t263, 0, -t261, 0; 0, 0, -t265 * t280 - t277, 0, -t263 * t279 - t265 * t282, 0; t261, 0, t269 * t275 + (-t264 * t286 + t266 * t279) * t263, 0, t260, 0; t259, 0, t269 * t276 + (t264 * t279 + t266 * t286) * t263, 0, t262, 0; 0, 0, -t265 * t279 + t278, 0, t263 * t280 - t265 * t281, 0; -t263 * t287 - t276, 0, -t263 * t283 - t264 * t288, 0, 0, 0; -qJD(1) * t289 + t275, 0, -qJD(3) * t289 + t265 * t287, 0, 0, 0; 0, 0, t284, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:48:00
	% EndTime: 2019-10-10 00:48:00
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (333->32), mult. (233->56), div. (0->0), fcn. (233->6), ass. (0->39)
	t329 = qJ(3) + pkin(11);
	t322 = sin(t329);
	t330 = qJ(1) + pkin(10);
	t323 = sin(t330);
	t353 = t322 * t323;
	t331 = qJ(5) + qJ(6);
	t326 = sin(t331);
	t328 = qJD(5) + qJD(6);
	t352 = t326 * t328;
	t327 = cos(t331);
	t351 = t327 * t328;
	t324 = cos(t329);
	t350 = qJD(1) * t324;
	t325 = cos(t330);
	t349 = qJD(1) * t325;
	t348 = qJD(1) * t326;
	t347 = qJD(1) * t327;
	t346 = qJD(3) * t324;
	t345 = qJD(3) * t325;
	t344 = qJD(3) * t326;
	t343 = qJD(3) * t327;
	t342 = t322 * t344;
	t341 = t322 * t343;
	t340 = qJD(3) * t353;
	t339 = t323 * t346;
	t338 = t324 * t345;
	t337 = t324 * t328 - qJD(1);
	t336 = -t328 + t350;
	t335 = t323 * t336;
	t334 = t323 * t351 + t325 * t348;
	t333 = t323 * t347 + t325 * t352;
	t332 = t337 * t326 + t341;
	t318 = t322 * t352 - t324 * t343;
	t317 = -t322 * t351 - t324 * t344;
	t316 = -t336 * t327 * t325 + t332 * t323;
	t315 = t334 * t324 - t326 * t340 - t333;
	t314 = t332 * t325 + t327 * t335;
	t313 = t326 * t335 + (-t337 * t327 + t342) * t325;
	t1 = [t316, 0, t333 * t322 - t327 * t338, 0, t313, t313; -t314, 0, -t327 * t339 + (t323 * t352 - t325 * t347) * t322, 0, -t315, -t315; 0, 0, -t324 * t352 - t341, 0, t317, t317; t315, 0, t326 * t338 + (-t323 * t348 + t325 * t351) * t322, 0, t314, t314; t313, 0, t334 * t322 + t326 * t339, 0, t316, t316; 0, 0, -t324 * t351 + t342, 0, t318, t318; -t322 * t349 - t339, 0, -t322 * t345 - t323 * t350, 0, 0, 0; -qJD(1) * t353 + t338, 0, t324 * t349 - t340, 0, 0, 0; 0, 0, t346, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end