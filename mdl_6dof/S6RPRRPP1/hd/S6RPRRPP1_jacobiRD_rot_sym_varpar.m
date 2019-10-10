% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
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
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(9);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 01:09:43
	% EndTime: 2019-10-10 01:09:43
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
	t242 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 01:09:43
	% EndTime: 2019-10-10 01:09:43
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (168->27), mult. (173->43), div. (0->0), fcn. (173->6), ass. (0->33)
	t265 = qJ(4) + pkin(10);
	t263 = cos(t265);
	t266 = qJ(1) + pkin(9);
	t264 = cos(t266);
	t288 = t263 * t264;
	t267 = sin(qJ(3));
	t287 = qJD(1) * t267;
	t268 = cos(qJ(3));
	t286 = qJD(1) * t268;
	t285 = qJD(3) * t267;
	t284 = qJD(3) * t268;
	t283 = qJD(4) * t267;
	t282 = qJD(4) * t268;
	t281 = t264 * t287;
	t280 = t263 * t285;
	t262 = sin(t266);
	t279 = t262 * t285;
	t278 = t264 * t285;
	t261 = sin(t265);
	t277 = t261 * t283;
	t276 = t263 * t283;
	t275 = -qJD(1) + t282;
	t274 = -qJD(4) + t286;
	t273 = t275 * t261;
	t272 = t262 * t284 + t281;
	t271 = t262 * t287 - t264 * t284;
	t270 = -t263 * t284 + t277;
	t269 = t274 * t262 + t278;
	t260 = -t274 * t288 + (t273 + t280) * t262;
	t259 = t275 * t263 * t262 + (t274 * t264 - t279) * t261;
	t258 = t269 * t263 + t264 * t273;
	t257 = t269 * t261 - t275 * t288;
	t1 = [t260, 0, t271 * t263 + t264 * t277, t257, 0, 0; -t258, 0, t270 * t262 - t263 * t281, -t259, 0, 0; 0, 0, -t261 * t282 - t280, -t261 * t284 - t276, 0, 0; t259, 0, -t271 * t261 + t264 * t276, t258, 0, 0; t257, 0, t272 * t261 + t262 * t276, t260, 0, 0; 0, 0, t261 * t285 - t263 * t282, t270, 0, 0; -t272, 0, -t262 * t286 - t278, 0, 0, 0; -t271, 0, t264 * t286 - t279, 0, 0, 0; 0, 0, t284, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:44
	% EndTime: 2019-10-10 01:09:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (168->26), mult. (173->43), div. (0->0), fcn. (173->6), ass. (0->33)
	t329 = qJ(1) + pkin(9);
	t327 = cos(t329);
	t331 = cos(qJ(3));
	t348 = qJD(1) * t331;
	t336 = -qJD(4) + t348;
	t351 = t336 * t327;
	t325 = sin(t329);
	t330 = sin(qJ(3));
	t347 = qJD(3) * t330;
	t340 = t327 * t347;
	t350 = t336 * t325 + t340;
	t349 = qJD(1) * t330;
	t346 = qJD(3) * t331;
	t345 = qJD(4) * t330;
	t344 = qJD(4) * t331;
	t343 = t327 * t349;
	t328 = qJ(4) + pkin(10);
	t326 = cos(t328);
	t342 = t326 * t347;
	t341 = t325 * t347;
	t324 = sin(t328);
	t339 = t324 * t345;
	t338 = t326 * t345;
	t337 = qJD(1) - t344;
	t335 = t337 * t327;
	t334 = -t325 * t346 - t343;
	t333 = t325 * t349 - t327 * t346;
	t332 = t326 * t346 - t339;
	t323 = t326 * t351 + (t337 * t324 - t342) * t325;
	t322 = t337 * t326 * t325 + (t341 - t351) * t324;
	t321 = t324 * t335 - t350 * t326;
	t320 = t350 * t324 + t326 * t335;
	t1 = [-t323, 0, t333 * t326 + t327 * t339, t320, 0, 0; t321, 0, -t332 * t325 - t326 * t343, t322, 0, 0; 0, 0, -t324 * t344 - t342, -t324 * t346 - t338, 0, 0; t334, 0, -t325 * t348 - t340, 0, 0, 0; -t333, 0, t327 * t348 - t341, 0, 0, 0; 0, 0, t346, 0, 0, 0; t322, 0, t333 * t324 - t327 * t338, t321, 0, 0; -t320, 0, t334 * t324 - t325 * t338, t323, 0, 0; 0, 0, -t324 * t347 + t326 * t344, t332, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end