% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (13->13), mult. (43->29), div. (0->0), fcn. (47->11), ass. (0->15)
	t159 = sin(pkin(13));
	t166 = cos(pkin(6));
	t172 = t159 * t166;
	t161 = sin(pkin(7));
	t162 = sin(pkin(6));
	t171 = t162 * t161;
	t164 = cos(pkin(13));
	t170 = t164 * t166;
	t169 = qJD(3) * sin(pkin(8));
	t168 = cos(qJ(3));
	t167 = sin(qJ(3));
	t165 = cos(pkin(7));
	t163 = cos(pkin(14));
	t158 = sin(pkin(14));
	t1 = [0, 0, 0, -(-(-t158 * t172 + t164 * t163) * t168 + (-(-t164 * t158 - t163 * t172) * t165 - t159 * t171) * t167) * t169, 0, 0; 0, 0, 0, -(-(t158 * t170 + t159 * t163) * t168 + (-(-t159 * t158 + t163 * t170) * t165 + t164 * t171) * t167) * t169, 0, 0; 0, 0, 0, -(-t161 * t166 * t167 + (-t163 * t165 * t167 - t158 * t168) * t162) * t169, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:22
	% EndTime: 2019-10-09 21:22:22
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (70->30), mult. (233->69), div. (0->0), fcn. (270->14), ass. (0->38)
	t242 = sin(pkin(14));
	t246 = sin(pkin(6));
	t253 = sin(qJ(3));
	t255 = cos(qJ(3));
	t247 = cos(pkin(14));
	t250 = cos(pkin(7));
	t266 = t247 * t250;
	t245 = sin(pkin(7));
	t251 = cos(pkin(6));
	t268 = t245 * t251;
	t275 = (t242 * t255 + t253 * t266) * t246 + t253 * t268;
	t248 = cos(pkin(13));
	t243 = sin(pkin(13));
	t270 = t243 * t251;
	t241 = -t242 * t270 + t248 * t247;
	t240 = -t248 * t242 - t247 * t270;
	t269 = t245 * t246;
	t260 = t240 * t250 + t243 * t269;
	t274 = t241 * t255 + t260 * t253;
	t265 = t248 * t251;
	t239 = t242 * t265 + t243 * t247;
	t238 = -t243 * t242 + t247 * t265;
	t261 = -t238 * t250 + t248 * t269;
	t273 = -t239 * t255 + t261 * t253;
	t267 = t246 * t250;
	t249 = cos(pkin(8));
	t254 = cos(qJ(4));
	t264 = t249 * t254;
	t252 = sin(qJ(4));
	t263 = qJD(3) * t252;
	t258 = -t239 * t253 - t261 * t255;
	t257 = -t241 * t253 + t260 * t255;
	t256 = t255 * t268 + (-t242 * t253 + t255 * t266) * t246;
	t244 = sin(pkin(8));
	t237 = t275 * qJD(3);
	t236 = t274 * qJD(3);
	t235 = t273 * qJD(3);
	t1 = [0, 0, 0, t236 * t244, t236 * t264 + (t274 * t254 + (t257 * t249 + (-t240 * t245 + t243 * t267) * t244) * t252) * qJD(4) + t257 * t263, 0; 0, 0, 0, -t235 * t244, -t235 * t264 + (-t273 * t254 + (t258 * t249 + (-t238 * t245 - t248 * t267) * t244) * t252) * qJD(4) + t258 * t263, 0; 0, 0, 0, t237 * t244, t237 * t264 + (t275 * t254 + (t256 * t249 + (-t247 * t269 + t251 * t250) * t244) * t252) * qJD(4) + t256 * t263, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:22
	% EndTime: 2019-10-09 21:22:23
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (195->51), mult. (628->110), div. (0->0), fcn. (751->16), ass. (0->54)
	t330 = sin(pkin(14));
	t334 = sin(pkin(6));
	t342 = sin(qJ(3));
	t345 = cos(qJ(3));
	t335 = cos(pkin(14));
	t338 = cos(pkin(7));
	t360 = t335 * t338;
	t333 = sin(pkin(7));
	t339 = cos(pkin(6));
	t362 = t333 * t339;
	t322 = (t330 * t345 + t342 * t360) * t334 + t342 * t362;
	t336 = cos(pkin(13));
	t331 = sin(pkin(13));
	t365 = t331 * t339;
	t329 = -t330 * t365 + t335 * t336;
	t328 = -t330 * t336 - t335 * t365;
	t363 = t333 * t334;
	t350 = t328 * t338 + t331 * t363;
	t318 = t329 * t345 + t350 * t342;
	t359 = t336 * t339;
	t327 = t330 * t359 + t331 * t335;
	t326 = -t330 * t331 + t335 * t359;
	t351 = -t326 * t338 + t336 * t363;
	t368 = -t327 * t345 + t351 * t342;
	t332 = sin(pkin(8));
	t343 = cos(qJ(5));
	t364 = t332 * t343;
	t361 = t334 * t338;
	t337 = cos(pkin(8));
	t341 = sin(qJ(4));
	t358 = t337 * t341;
	t344 = cos(qJ(4));
	t357 = t337 * t344;
	t340 = sin(qJ(5));
	t356 = qJD(4) * t340;
	t315 = -t327 * t342 - t351 * t345;
	t323 = -t326 * t333 - t336 * t361;
	t354 = t315 * t337 + t323 * t332;
	t317 = -t329 * t342 + t350 * t345;
	t324 = -t328 * t333 + t331 * t361;
	t353 = t317 * t337 + t324 * t332;
	t321 = t345 * t362 + (-t330 * t342 + t345 * t360) * t334;
	t325 = -t335 * t363 + t338 * t339;
	t352 = t321 * t337 + t325 * t332;
	t348 = t354 * t341 - t344 * t368;
	t347 = t318 * t344 + t353 * t341;
	t346 = t322 * t344 + t352 * t341;
	t320 = t322 * qJD(3);
	t319 = t321 * qJD(3);
	t314 = t318 * qJD(3);
	t313 = t317 * qJD(3);
	t312 = t368 * qJD(3);
	t311 = t315 * qJD(3);
	t1 = [0, 0, 0, t314 * t332, t347 * qJD(4) + t313 * t341 + t314 * t357, (t313 * t344 - t314 * t358) * t340 - t314 * t364 + (t347 * t343 + (-t317 * t332 + t324 * t337) * t340) * qJD(5) + (-t318 * t341 + t353 * t344) * t356; 0, 0, 0, -t312 * t332, t348 * qJD(4) + t311 * t341 - t312 * t357, (t311 * t344 + t312 * t358) * t340 + t312 * t364 + (t348 * t343 + (-t315 * t332 + t323 * t337) * t340) * qJD(5) + (t341 * t368 + t354 * t344) * t356; 0, 0, 0, t320 * t332, t346 * qJD(4) + t319 * t341 + t320 * t357, (t319 * t344 - t320 * t358) * t340 - t320 * t364 + (t346 * t343 + (-t321 * t332 + t325 * t337) * t340) * qJD(5) + (-t322 * t341 + t352 * t344) * t356;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end