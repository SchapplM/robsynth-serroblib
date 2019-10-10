% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(11)) + r_i_i_C(2) * sin(pkin(11)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(11) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(11)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (124->30), mult. (120->40), div. (0->0), fcn. (79->7), ass. (0->28)
	t44 = qJD(3) + qJD(4);
	t43 = pkin(11) + qJ(3);
	t41 = qJ(4) + t43;
	t38 = cos(t41);
	t61 = r_i_i_C(2) * t38;
	t37 = sin(t41);
	t63 = r_i_i_C(1) * t37;
	t51 = t61 + t63;
	t49 = t51 * t44;
	t39 = sin(t43);
	t64 = pkin(3) * t39;
	t65 = qJD(3) * t64 + t49;
	t62 = r_i_i_C(2) * t37;
	t60 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
	t59 = t38 * t44;
	t45 = sin(qJ(1));
	t58 = qJD(1) * t45;
	t46 = cos(qJ(1));
	t57 = qJD(1) * t46;
	t40 = cos(t43);
	t56 = qJD(3) * t40;
	t55 = r_i_i_C(1) * t59;
	t54 = t44 * t62;
	t52 = qJD(1) * t61;
	t50 = -r_i_i_C(1) * t38 - pkin(3) * t40 - cos(pkin(11)) * pkin(2) - pkin(1) + t62;
	t48 = t45 * t52 + t58 * t63 + (t54 - t55) * t46;
	t32 = t45 * t54;
	t1 = [t46 * qJD(2) + t65 * t45 + (-t60 * t45 + t50 * t46) * qJD(1), t57, (t39 * t58 - t46 * t56) * pkin(3) + t48, t48, 0, 0; t45 * qJD(2) - t65 * t46 + (t50 * t45 + t60 * t46) * qJD(1), t58, t32 + (-pkin(3) * t56 - t55) * t45 + (-t51 - t64) * t57, -t46 * t52 + t32 + (-t37 * t57 - t45 * t59) * r_i_i_C(1), 0, 0; 0, 0, -t65, -t49, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:08
	% EndTime: 2019-10-10 09:06:09
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (387->62), mult. (400->98), div. (0->0), fcn. (303->9), ass. (0->55)
	t262 = pkin(11) + qJ(3);
	t260 = qJ(4) + t262;
	t256 = sin(t260);
	t266 = cos(qJ(5));
	t312 = r_i_i_C(1) * t266 + pkin(4);
	t315 = t256 * t312;
	t264 = sin(qJ(5));
	t295 = qJD(5) * t266;
	t257 = cos(t260);
	t263 = qJD(3) + qJD(4);
	t303 = t257 * t263;
	t314 = t256 * t295 + t264 * t303;
	t309 = pkin(9) + r_i_i_C(3);
	t290 = t309 * t257;
	t258 = sin(t262);
	t304 = pkin(3) * qJD(3);
	t293 = t258 * t304;
	t307 = pkin(4) * t256;
	t313 = -t293 + (t290 - t307) * t263;
	t296 = qJD(5) * t264;
	t283 = t256 * t296;
	t310 = r_i_i_C(1) * t283 + r_i_i_C(2) * t314;
	t308 = pkin(3) * t258;
	t305 = r_i_i_C(2) * t264;
	t265 = sin(qJ(1));
	t302 = t263 * t265;
	t301 = t263 * t266;
	t267 = cos(qJ(1));
	t300 = t263 * t267;
	t299 = t266 * t267;
	t298 = qJD(1) * t265;
	t297 = qJD(1) * t267;
	t294 = t256 * t305;
	t292 = qJD(1) * t305;
	t291 = t309 * t256;
	t289 = t309 * t265;
	t288 = t256 * t301;
	t278 = qJD(5) * t257 - qJD(1);
	t277 = qJD(1) * t257 - qJD(5);
	t276 = t312 * t257;
	t275 = t312 * t267;
	t274 = t267 * t310 + t298 * t315;
	t273 = t278 * t264;
	t272 = t267 * t256 * t292 + t265 * t310 + t290 * t297;
	t259 = cos(t262);
	t271 = -pkin(4) * t257 - pkin(3) * t259 - cos(pkin(11)) * pkin(2) - pkin(1) - t291;
	t270 = t256 * t300 + t265 * t277;
	t269 = -t259 * t304 + (-t276 - t291) * t263;
	t268 = -t257 * r_i_i_C(2) * t295 + (-t257 * t296 - t288) * r_i_i_C(1) + t309 * t303 + (-t307 + t294) * t263;
	t261 = -pkin(8) - pkin(7) - qJ(2);
	t240 = -t277 * t299 + (t273 + t288) * t265;
	t239 = t278 * t266 * t265 + (-t256 * t302 + t267 * t277) * t264;
	t238 = t266 * t270 + t267 * t273;
	t237 = t264 * t270 - t278 * t299;
	t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t267 * qJD(2) - t313 * t265 + (t261 * t265 + t267 * t271) * qJD(1), t297, (-t290 - t294 + t308) * t298 + t269 * t267 + t274, (-t265 * t292 - t300 * t309) * t256 + (-qJD(1) * t289 - t263 * t275) * t257 + t274, r_i_i_C(1) * t237 + r_i_i_C(2) * t238, 0; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t265 * qJD(2) + t313 * t267 + (-t261 * t267 + t265 * t271) * qJD(1), t298, (-t308 - t315) * t297 + t269 * t265 + t272, -t276 * t302 + (-qJD(1) * t275 - t263 * t289) * t256 + t272, -r_i_i_C(1) * t239 + r_i_i_C(2) * t240, 0; 0, 0, t268 - t293, t268, (-t257 * t301 + t283) * r_i_i_C(2) - t314 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:08
	% EndTime: 2019-10-10 09:06:09
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (673->81), mult. (566->111), div. (0->0), fcn. (437->11), ass. (0->70)
	t295 = pkin(11) + qJ(3);
	t291 = qJ(4) + t295;
	t286 = sin(t291);
	t298 = qJ(5) + qJ(6);
	t292 = sin(t298);
	t293 = cos(t298);
	t296 = qJD(5) + qJD(6);
	t345 = t293 * t296;
	t287 = cos(t291);
	t297 = qJD(3) + qJD(4);
	t347 = t287 * t297;
	t367 = t286 * t345 + t292 * t347;
	t303 = -pkin(10) - pkin(9);
	t299 = sin(qJ(5));
	t350 = pkin(5) * qJD(5);
	t335 = t299 * t350;
	t366 = t297 * t303 + t335;
	t301 = cos(qJ(5));
	t288 = t301 * pkin(5) + pkin(4);
	t364 = r_i_i_C(1) * t293 + t288;
	t365 = t286 * t364 + t287 * t303;
	t289 = sin(t295);
	t351 = pkin(3) * qJD(3);
	t336 = t289 * t351;
	t348 = t286 * t297;
	t352 = r_i_i_C(3) - t303;
	t356 = pkin(5) * t299;
	t363 = (t352 * t297 - t335) * t287 + (pkin(8) + pkin(7) + qJ(2) + t356) * qJD(1) - t288 * t348 - t336;
	t349 = t286 * t292;
	t332 = t296 * t349;
	t362 = r_i_i_C(1) * t332 + t367 * r_i_i_C(2) + t366 * t286;
	t302 = cos(qJ(1));
	t320 = t287 * t296 - qJD(1);
	t361 = t302 * t320;
	t340 = qJD(1) * t287;
	t319 = -t296 + t340;
	t300 = sin(qJ(1));
	t330 = t300 * t348;
	t359 = t319 * t302 - t330;
	t357 = pkin(3) * t289;
	t354 = r_i_i_C(2) * t293;
	t353 = r_i_i_C(3) * t287;
	t344 = t297 * t302;
	t329 = t286 * t344;
	t307 = t319 * t300 + t329;
	t263 = t307 * t292 - t293 * t361;
	t264 = t292 * t361 + t307 * t293;
	t342 = t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
	t314 = t320 * t300;
	t265 = t359 * t292 + t293 * t314;
	t266 = t292 * t314 - t359 * t293;
	t341 = -t265 * r_i_i_C(1) + t266 * r_i_i_C(2);
	t339 = qJD(1) * t300;
	t338 = qJD(1) * t302;
	t337 = r_i_i_C(2) * t349;
	t334 = t301 * t350;
	t333 = r_i_i_C(2) * qJD(1) * t292;
	t317 = -qJD(5) + t340;
	t316 = -r_i_i_C(1) * t292 - t354;
	t315 = t364 * t297;
	t313 = (-qJD(5) * t287 + qJD(1)) * t301;
	t312 = t302 * t286 * t333 + t362 * t300 + t338 * t353;
	t309 = t362 * t302 + t365 * t339;
	t308 = (-r_i_i_C(3) * t286 - t287 * t364) * t297;
	t290 = cos(t295);
	t306 = -t290 * t351 + t308;
	t305 = t334 + qJD(2) + (-t352 * t286 - t287 * t288 - pkin(3) * t290 - cos(pkin(11)) * pkin(2) - pkin(1)) * qJD(1);
	t304 = t297 * t337 + r_i_i_C(3) * t347 - t286 * t315 + (t316 * t296 - t366) * t287;
	t276 = r_i_i_C(2) * t332;
	t1 = [t266 * r_i_i_C(1) + t265 * r_i_i_C(2) - t363 * t300 + t305 * t302, t338, (-t337 - t353 + t357) * t339 + t306 * t302 + t309, (-r_i_i_C(3) * t344 - t300 * t333) * t286 + (-r_i_i_C(3) * t339 - t302 * t315) * t287 + t309, (t302 * t313 + (t317 * t300 + t329) * t299) * pkin(5) + t342, t342; -t264 * r_i_i_C(1) + t263 * r_i_i_C(2) + t305 * t300 + t363 * t302, t339, t306 * t300 + (-t365 - t357) * t338 + t312, t300 * t308 - t338 * t365 + t312, (t300 * t313 + (-t317 * t302 + t330) * t299) * pkin(5) + t341, t341; 0, 0, t304 - t336, t304, t276 + (-r_i_i_C(1) * t345 - t334) * t286 + (t316 - t356) * t347, -t367 * r_i_i_C(1) - t347 * t354 + t276;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end