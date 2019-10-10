% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:30
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
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
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(10);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(7);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (131->30), mult. (132->40), div. (0->0), fcn. (86->8), ass. (0->27)
	t49 = qJ(2) + pkin(10);
	t41 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t49);
	t48 = qJD(2) + qJD(4);
	t46 = qJ(4) + t49;
	t43 = cos(t46);
	t67 = r_i_i_C(2) * t43;
	t42 = sin(t46);
	t69 = r_i_i_C(1) * t42;
	t56 = t67 + t69;
	t54 = t56 * t48;
	t70 = t41 * qJD(2) - t54;
	t68 = r_i_i_C(2) * t42;
	t66 = r_i_i_C(3) + pkin(8) + qJ(3) + pkin(7);
	t65 = t43 * t48;
	t51 = sin(qJ(1));
	t64 = qJD(1) * t51;
	t53 = cos(qJ(1));
	t63 = qJD(1) * t53;
	t62 = r_i_i_C(1) * t65;
	t61 = t48 * t68;
	t59 = qJD(1) * t67;
	t60 = t51 * t59 + t53 * t61 + t64 * t69;
	t57 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t49);
	t58 = t57 * qJD(2) - t62;
	t55 = -r_i_i_C(1) * t43 - pkin(1) + t57 + t68;
	t36 = t51 * t61;
	t1 = [t53 * qJD(3) - t70 * t51 + (-t66 * t51 + t55 * t53) * qJD(1), -t41 * t64 + t58 * t53 + t60, t63, -t53 * t62 + t60, 0, 0; t51 * qJD(3) + t70 * t53 + (t55 * t51 + t66 * t53) * qJD(1), t36 + t58 * t51 + (t41 - t56) * t63, t64, -t53 * t59 + t36 + (-t42 * t63 - t51 * t65) * r_i_i_C(1), 0, 0; 0, t70, 0, -t54, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:07
	% EndTime: 2019-10-10 10:30:07
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (394->63), mult. (412->98), div. (0->0), fcn. (310->10), ass. (0->53)
	t268 = qJ(2) + pkin(10);
	t265 = qJ(4) + t268;
	t261 = sin(t265);
	t272 = cos(qJ(5));
	t317 = r_i_i_C(1) * t272 + pkin(4);
	t320 = t261 * t317;
	t269 = sin(qJ(5));
	t302 = qJD(5) * t272;
	t262 = cos(t265);
	t267 = qJD(2) + qJD(4);
	t310 = t262 * t267;
	t319 = t261 * t302 + t269 * t310;
	t314 = pkin(9) + r_i_i_C(3);
	t298 = t314 * t262;
	t258 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t268);
	t251 = t258 * qJD(2);
	t313 = pkin(4) * t261;
	t318 = t251 + (t298 - t313) * t267;
	t303 = qJD(5) * t269;
	t291 = t261 * t303;
	t315 = r_i_i_C(1) * t291 + r_i_i_C(2) * t319;
	t311 = r_i_i_C(2) * t269;
	t271 = sin(qJ(1));
	t309 = t267 * t271;
	t308 = t267 * t272;
	t274 = cos(qJ(1));
	t307 = t267 * t274;
	t306 = t272 * t274;
	t305 = qJD(1) * t271;
	t304 = qJD(1) * t274;
	t301 = t261 * t311;
	t300 = qJD(1) * t311;
	t299 = t314 * t261;
	t297 = t314 * t271;
	t296 = t261 * t308;
	t286 = qJD(5) * t262 - qJD(1);
	t285 = qJD(1) * t262 - qJD(5);
	t284 = t317 * t262;
	t283 = t317 * t274;
	t282 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t268);
	t281 = t274 * t315 + t305 * t320;
	t280 = t286 * t269;
	t279 = t274 * t261 * t300 + t271 * t315 + t298 * t304;
	t278 = -pkin(4) * t262 - pkin(1) + t282 - t299;
	t277 = t261 * t307 + t271 * t285;
	t276 = t282 * qJD(2) + (-t284 - t299) * t267;
	t275 = -t262 * r_i_i_C(2) * t302 + (-t262 * t303 - t296) * r_i_i_C(1) + t314 * t310 + (-t313 + t301) * t267;
	t266 = -pkin(8) - qJ(3) - pkin(7);
	t242 = -t285 * t306 + (t280 + t296) * t271;
	t241 = t286 * t272 * t271 + (-t261 * t309 + t274 * t285) * t269;
	t240 = t272 * t277 + t274 * t280;
	t239 = t269 * t277 - t286 * t306;
	t1 = [t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t274 * qJD(3) - t318 * t271 + (t266 * t271 + t274 * t278) * qJD(1), (-t258 - t298 - t301) * t305 + t276 * t274 + t281, t304, (-t271 * t300 - t307 * t314) * t261 + (-qJD(1) * t297 - t267 * t283) * t262 + t281, r_i_i_C(1) * t239 + r_i_i_C(2) * t240, 0; -t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t271 * qJD(3) + t318 * t274 + (-t266 * t274 + t271 * t278) * qJD(1), (t258 - t320) * t304 + t276 * t271 + t279, t305, -t284 * t309 + (-qJD(1) * t283 - t267 * t297) * t261 + t279, -r_i_i_C(1) * t241 + r_i_i_C(2) * t242, 0; 0, t251 + t275, 0, t275, (-t262 * t308 + t291) * r_i_i_C(2) - t319 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:07
	% EndTime: 2019-10-10 10:30:07
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (520->75), mult. (533->105), div. (0->0), fcn. (407->10), ass. (0->55)
	t329 = pkin(5) + r_i_i_C(1);
	t272 = qJD(2) + qJD(4);
	t275 = sin(qJ(5));
	t323 = r_i_i_C(2) * t275;
	t336 = t272 * t323 + qJD(6);
	t273 = qJ(2) + pkin(10);
	t270 = qJ(4) + t273;
	t265 = sin(t270);
	t266 = cos(t270);
	t278 = cos(qJ(5));
	t312 = qJD(5) * t278;
	t313 = qJD(5) * t275;
	t335 = t336 * t266 + (r_i_i_C(2) * t312 + t313 * t329) * t265;
	t277 = sin(qJ(1));
	t280 = cos(qJ(1));
	t297 = qJD(1) * t266 - qJD(5);
	t290 = t297 * t280;
	t298 = qJD(5) * t266 - qJD(1);
	t292 = t298 * t278;
	t317 = t272 * t277;
	t308 = t265 * t317;
	t239 = (t290 - t308) * t275 + t277 * t292;
	t267 = t278 * pkin(5) + pkin(4);
	t334 = r_i_i_C(1) * t278 + t267;
	t260 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t273);
	t254 = t260 * qJD(2);
	t274 = -qJ(6) - pkin(9);
	t321 = r_i_i_C(3) - t274;
	t330 = (-pkin(5) * t313 + t321 * t272) * t266 + (pkin(5) * t275 + pkin(7) + pkin(8) + qJ(3)) * qJD(1) - (t267 * t272 - qJD(6)) * t265 + t254;
	t322 = r_i_i_C(3) * t266;
	t320 = t266 * t272;
	t319 = t266 * t274;
	t318 = t272 * t265;
	t316 = t272 * t280;
	t315 = qJD(1) * t277;
	t314 = qJD(1) * t280;
	t307 = t265 * t316;
	t306 = t265 * t315;
	t305 = t265 * t314;
	t294 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t273);
	t293 = t334 * t272;
	t291 = t298 * t275;
	t289 = -r_i_i_C(2) * t278 - t275 * t329;
	t288 = -t265 * t334 - t319;
	t287 = t274 * t308 + t335 * t277 + t305 * t323 + t314 * t322;
	t286 = t274 * t307 + t335 * t280 + t334 * t306 + t315 * t319;
	t285 = (-r_i_i_C(3) * t265 - t266 * t334) * t272;
	t284 = t297 * t277 + t307;
	t283 = t294 * qJD(2) + t285;
	t282 = pkin(5) * t312 + qJD(3) + (-t321 * t265 - t266 * t267 - pkin(1) + t294) * qJD(1);
	t237 = t284 * t275 - t280 * t292;
	t281 = r_i_i_C(3) * t320 + (t289 * qJD(5) - t272 * t274) * t266 + (-t293 + t336) * t265;
	t240 = -t278 * t290 + (t278 * t318 + t291) * t277;
	t238 = t284 * t278 + t280 * t291;
	t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t330 * t277 + t282 * t280, (-t265 * t323 - t260 - t322) * t315 + t283 * t280 + t286, t314, (-r_i_i_C(3) * t316 - t315 * t323) * t265 + (-r_i_i_C(3) * t315 - t280 * t293) * t266 + t286, t238 * r_i_i_C(2) + t329 * t237, t266 * t316 - t306; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t282 * t277 + t330 * t280, t283 * t277 + (t260 + t288) * t314 + t287, t315, t277 * t285 + t288 * t314 + t287, t240 * r_i_i_C(2) - t329 * t239, t266 * t317 + t305; 0, t254 + t281, 0, t281, t289 * t320 + (-t278 * t329 + t323) * t265 * qJD(5), t318;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end