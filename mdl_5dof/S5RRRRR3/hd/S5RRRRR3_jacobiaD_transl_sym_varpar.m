% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (15->13), mult. (56->29), div. (0->0), fcn. (36->4), ass. (0->12)
	t16 = sin(qJ(1));
	t25 = qJD(1) * t16;
	t18 = cos(qJ(1));
	t24 = qJD(1) * t18;
	t23 = qJD(2) * t16;
	t22 = qJD(2) * t18;
	t15 = sin(qJ(2));
	t17 = cos(qJ(2));
	t21 = -r_i_i_C(1) * t17 + r_i_i_C(2) * t15;
	t20 = r_i_i_C(1) * t15 + r_i_i_C(2) * t17;
	t19 = t20 * qJD(2);
	t1 = [t20 * t23 + (-r_i_i_C(3) * t16 + t21 * t18) * qJD(1), (t15 * t22 + t17 * t25) * r_i_i_C(2) + (t15 * t25 - t17 * t22) * r_i_i_C(1), 0, 0, 0; -t18 * t19 + (r_i_i_C(3) * t18 + t21 * t16) * qJD(1), (t15 * t23 - t17 * t24) * r_i_i_C(2) + (-t15 * t24 - t17 * t23) * r_i_i_C(1), 0, 0, 0; 0, -t19, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (75->24), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t34 = qJD(2) + qJD(3);
	t35 = qJ(2) + qJ(3);
	t33 = cos(t35);
	t53 = r_i_i_C(2) * t33;
	t32 = sin(t35);
	t55 = r_i_i_C(1) * t32;
	t44 = t53 + t55;
	t43 = t44 * t34;
	t36 = sin(qJ(2));
	t56 = pkin(1) * t36;
	t57 = qJD(2) * t56 + t43;
	t54 = r_i_i_C(2) * t32;
	t52 = t33 * t34;
	t37 = sin(qJ(1));
	t51 = qJD(1) * t37;
	t39 = cos(qJ(1));
	t50 = qJD(1) * t39;
	t38 = cos(qJ(2));
	t49 = qJD(2) * t38;
	t48 = r_i_i_C(1) * t52;
	t47 = t34 * t54;
	t46 = qJD(1) * t53;
	t42 = -pkin(1) * t38 - r_i_i_C(1) * t33 + t54;
	t41 = t37 * t46 + t51 * t55 + (t47 - t48) * t39;
	t28 = t37 * t47;
	t1 = [t57 * t37 + (-r_i_i_C(3) * t37 + t42 * t39) * qJD(1), (t36 * t51 - t39 * t49) * pkin(1) + t41, t41, 0, 0; -t57 * t39 + (r_i_i_C(3) * t39 + t42 * t37) * qJD(1), t28 + (-pkin(1) * t49 - t48) * t37 + (-t44 - t56) * t50, -t39 * t46 + t28 + (-t32 * t50 - t37 * t52) * r_i_i_C(1), 0, 0; 0, -t57, -t43, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:42
	% EndTime: 2019-10-09 21:04:42
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (264->54), mult. (390->92), div. (0->0), fcn. (295->8), ass. (0->53)
	t256 = qJ(2) + qJ(3);
	t253 = sin(t256);
	t260 = cos(qJ(4));
	t307 = r_i_i_C(1) * t260 + pkin(2);
	t310 = t253 * t307;
	t257 = sin(qJ(4));
	t290 = qJD(4) * t260;
	t254 = cos(t256);
	t255 = qJD(2) + qJD(3);
	t298 = t254 * t255;
	t309 = t253 * t290 + t257 * t298;
	t304 = pkin(5) + r_i_i_C(3);
	t285 = t304 * t254;
	t258 = sin(qJ(2));
	t299 = pkin(1) * qJD(2);
	t287 = t258 * t299;
	t302 = pkin(2) * t253;
	t308 = -t287 + (t285 - t302) * t255;
	t291 = qJD(4) * t257;
	t278 = t253 * t291;
	t305 = r_i_i_C(1) * t278 + t309 * r_i_i_C(2);
	t303 = pkin(1) * t258;
	t300 = r_i_i_C(2) * t257;
	t259 = sin(qJ(1));
	t297 = t255 * t259;
	t296 = t255 * t260;
	t262 = cos(qJ(1));
	t295 = t255 * t262;
	t294 = t260 * t262;
	t293 = qJD(1) * t259;
	t292 = qJD(1) * t262;
	t289 = t253 * t300;
	t288 = qJD(1) * t300;
	t286 = t304 * t253;
	t284 = t304 * t259;
	t283 = t253 * t296;
	t273 = qJD(4) * t254 - qJD(1);
	t272 = qJD(1) * t254 - qJD(4);
	t271 = t307 * t254;
	t270 = t307 * t262;
	t269 = t305 * t262 + t293 * t310;
	t268 = t273 * t257;
	t267 = t262 * t253 * t288 + t305 * t259 + t292 * t285;
	t266 = t253 * t295 + t272 * t259;
	t261 = cos(qJ(2));
	t265 = qJD(1) * (-pkin(1) * t261 - pkin(2) * t254 - t286);
	t264 = -t261 * t299 + (-t271 - t286) * t255;
	t263 = -t254 * r_i_i_C(2) * t290 + (-t254 * t291 - t283) * r_i_i_C(1) + t304 * t298 + (-t302 + t289) * t255;
	t238 = -t272 * t294 + (t268 + t283) * t259;
	t237 = t273 * t260 * t259 + (-t253 * t297 + t272 * t262) * t257;
	t236 = t266 * t260 + t262 * t268;
	t235 = t266 * t257 - t273 * t294;
	t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) - t308 * t259 + t262 * t265, (-t285 - t289 + t303) * t293 + t264 * t262 + t269, (-t259 * t288 - t304 * t295) * t253 + (-qJD(1) * t284 - t255 * t270) * t254 + t269, t235 * r_i_i_C(1) + t236 * r_i_i_C(2), 0; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t259 * t265 + t308 * t262, (-t303 - t310) * t292 + t264 * t259 + t267, -t271 * t297 + (-qJD(1) * t270 - t255 * t284) * t253 + t267, -t237 * r_i_i_C(1) + t238 * r_i_i_C(2), 0; 0, t263 - t287, t263, (-t254 * t296 + t278) * r_i_i_C(2) - t309 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:42
	% EndTime: 2019-10-09 21:04:42
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (504->72), mult. (556->111), div. (0->0), fcn. (429->10), ass. (0->70)
	t291 = qJ(2) + qJ(3);
	t285 = sin(t291);
	t295 = cos(qJ(4));
	t283 = t295 * pkin(3) + pkin(2);
	t290 = qJ(4) + qJ(5);
	t286 = cos(t290);
	t354 = r_i_i_C(1) * t286 + t283;
	t358 = t285 * t354;
	t284 = sin(t290);
	t287 = cos(t291);
	t289 = qJD(2) + qJD(3);
	t338 = t287 * t289;
	t288 = qJD(4) + qJD(5);
	t339 = t286 * t288;
	t357 = t284 * t338 + t285 * t339;
	t349 = pkin(5) + r_i_i_C(3);
	t356 = t349 * t287;
	t326 = t349 * t289;
	t293 = sin(qJ(2));
	t344 = pkin(1) * qJD(2);
	t328 = t293 * t344;
	t292 = sin(qJ(4));
	t343 = pkin(3) * qJD(4);
	t331 = t292 * t343;
	t341 = t285 * t289;
	t347 = pkin(3) * t292;
	t355 = qJD(1) * t347 - t283 * t341 + (t326 - t331) * t287 - t328;
	t342 = t284 * t285;
	t324 = t288 * t342;
	t352 = r_i_i_C(1) * t324 + t357 * r_i_i_C(2) + t285 * t331;
	t297 = cos(qJ(1));
	t313 = t287 * t288 - qJD(1);
	t351 = t297 * t313;
	t335 = qJD(1) * t287;
	t312 = -t288 + t335;
	t294 = sin(qJ(1));
	t321 = t294 * t341;
	t350 = t312 * t297 - t321;
	t348 = pkin(1) * t293;
	t345 = r_i_i_C(2) * t286;
	t340 = t285 * t297;
	t320 = t289 * t340;
	t302 = t312 * t294 + t320;
	t262 = t302 * t284 - t286 * t351;
	t263 = t284 * t351 + t302 * t286;
	t337 = t262 * r_i_i_C(1) + t263 * r_i_i_C(2);
	t306 = t313 * t294;
	t264 = t350 * t284 + t286 * t306;
	t265 = t284 * t306 - t350 * t286;
	t336 = -t264 * r_i_i_C(1) + t265 * r_i_i_C(2);
	t334 = qJD(1) * t294;
	t333 = qJD(1) * t297;
	t332 = r_i_i_C(2) * t342;
	t330 = t295 * t343;
	t329 = r_i_i_C(2) * qJD(1) * t284;
	t327 = t349 * t285;
	t325 = t349 * t294;
	t310 = -qJD(4) + t335;
	t309 = -r_i_i_C(1) * t284 - t345;
	t308 = t354 * t289;
	t307 = t354 * t297;
	t305 = (-qJD(4) * t287 + qJD(1)) * t295;
	t304 = t352 * t297 + t334 * t358;
	t303 = t352 * t294 + t329 * t340 + t333 * t356;
	t296 = cos(qJ(2));
	t300 = t330 + (-pkin(1) * t296 - t283 * t287 - t327) * qJD(1);
	t299 = -t296 * t344 + (-t287 * t354 - t327) * t289;
	t298 = t289 * t332 - t285 * t308 + (t309 * t288 - t331) * t287 + t349 * t338;
	t275 = r_i_i_C(2) * t324;
	t1 = [t265 * r_i_i_C(1) + t264 * r_i_i_C(2) - t355 * t294 + t300 * t297, (-t332 + t348 - t356) * t334 + t299 * t297 + t304, (-t294 * t329 - t297 * t326) * t285 + (-qJD(1) * t325 - t289 * t307) * t287 + t304, (t297 * t305 + (t310 * t294 + t320) * t292) * pkin(3) + t337, t337; -t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t300 * t294 + t355 * t297, (-t348 - t358) * t333 + t299 * t294 + t303, -t294 * t287 * t308 + (-qJD(1) * t307 - t289 * t325) * t285 + t303, (t294 * t305 + (-t310 * t297 + t321) * t292) * pkin(3) + t336, t336; 0, t298 - t328, t298, t275 + (-r_i_i_C(1) * t339 - t330) * t285 + (t309 - t347) * t338, -t357 * r_i_i_C(1) - t338 * t345 + t275;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end