% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
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
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
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
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:41
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:42
	% EndTime: 2019-10-10 11:38:42
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (217->43), mult. (285->60), div. (0->0), fcn. (208->8), ass. (0->38)
	t224 = qJ(2) + qJ(3);
	t222 = cos(t224);
	t259 = r_i_i_C(3) + qJ(4);
	t270 = t259 * t222;
	t221 = sin(t224);
	t219 = t221 * qJD(4);
	t223 = qJD(2) + qJD(3);
	t225 = sin(pkin(10));
	t226 = cos(pkin(10));
	t260 = r_i_i_C(2) * t225;
	t268 = r_i_i_C(1) * t226 + pkin(3);
	t237 = t268 - t260;
	t227 = sin(qJ(2));
	t258 = pkin(2) * qJD(2);
	t250 = t227 * t258;
	t269 = (-t237 * t221 + t270) * t223 + (t225 * r_i_i_C(1) + t226 * r_i_i_C(2) + pkin(7) + pkin(8)) * qJD(1) + t219 - t250;
	t251 = t223 * t260;
	t267 = (qJD(4) + t251) * t222;
	t262 = pkin(2) * t227;
	t228 = sin(qJ(1));
	t256 = t222 * t228;
	t230 = cos(qJ(1));
	t255 = t223 * t230;
	t254 = qJD(1) * t228;
	t253 = qJD(1) * t230;
	t248 = t221 * t254;
	t247 = t221 * t253;
	t245 = t259 * t221;
	t243 = t259 * t228;
	t241 = t267 * t230 + t268 * t248;
	t240 = t268 * t223;
	t239 = t268 * t230;
	t238 = t267 * t228 + t247 * t260 + t253 * t270;
	t234 = t219 + t223 * t270 + (-t240 + t251) * t221;
	t229 = cos(qJ(2));
	t233 = qJD(1) * (-t229 * pkin(2) - t237 * t222 - pkin(1) - t245);
	t232 = -t229 * t258 + (-t222 * t268 - t245) * t223;
	t1 = [-t269 * t228 + t230 * t233, (-t221 * t260 + t262 - t270) * t254 + t232 * t230 + t241, (-t254 * t260 - t259 * t255) * t221 + (-qJD(1) * t243 - t223 * t239) * t222 + t241, t222 * t255 - t248, 0, 0; t228 * t233 + t269 * t230, (-t221 * t268 - t262) * t253 + t232 * t228 + t238, -t240 * t256 + (-qJD(1) * t239 - t223 * t243) * t221 + t238, t223 * t256 + t247, 0, 0; 0, t234 - t250, t234, t223 * t221, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:42
	% EndTime: 2019-10-10 11:38:42
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (402->68), mult. (437->101), div. (0->0), fcn. (338->10), ass. (0->56)
	t259 = pkin(10) + qJ(5);
	t255 = sin(t259);
	t256 = cos(t259);
	t261 = qJ(2) + qJ(3);
	t257 = sin(t261);
	t296 = qJD(5) * t257;
	t258 = cos(t261);
	t260 = qJD(2) + qJD(3);
	t303 = t258 * t260;
	t318 = t255 * t303 + t256 * t296;
	t253 = cos(pkin(10)) * pkin(4) + pkin(3);
	t317 = r_i_i_C(1) * t256 + t253;
	t252 = t257 * qJD(4);
	t264 = sin(qJ(2));
	t305 = pkin(2) * qJD(2);
	t294 = t264 * t305;
	t263 = -pkin(9) - qJ(4);
	t306 = r_i_i_C(3) - t263;
	t316 = (-t253 * t257 + t306 * t258) * t260 + (pkin(4) * sin(pkin(10)) + pkin(8) + pkin(7)) * qJD(1) + t252 - t294;
	t288 = t255 * t296;
	t315 = r_i_i_C(1) * t288 + t318 * r_i_i_C(2) + qJD(4) * t258;
	t267 = cos(qJ(1));
	t281 = qJD(5) * t258 - qJD(1);
	t314 = t267 * t281;
	t280 = qJD(1) * t258 - qJD(5);
	t265 = sin(qJ(1));
	t301 = t260 * t265;
	t292 = t257 * t301;
	t312 = t280 * t267 - t292;
	t310 = pkin(2) * t264;
	t308 = r_i_i_C(2) * t255;
	t307 = r_i_i_C(3) * t258;
	t302 = t258 * t263;
	t300 = t260 * t267;
	t299 = qJD(1) * t265;
	t298 = qJD(1) * t267;
	t295 = t257 * t308;
	t291 = t257 * t300;
	t290 = t257 * t299;
	t289 = t257 * t298;
	t279 = t317 * t260;
	t278 = t281 * t265;
	t276 = t263 * t292 + t315 * t265 + t289 * t308 + t298 * t307;
	t275 = -t257 * t317 - t302;
	t274 = t263 * t291 + t315 * t267 + t317 * t290 + t299 * t302;
	t266 = cos(qJ(2));
	t273 = qJD(1) * (-t266 * pkin(2) - t253 * t258 - t306 * t257 - pkin(1));
	t272 = (-r_i_i_C(3) * t257 - t258 * t317) * t260;
	t271 = t280 * t265 + t291;
	t270 = -t266 * t305 + t272;
	t269 = t260 * t295 + r_i_i_C(3) * t303 + t252 - t257 * t279 + (-t260 * t263 + (-r_i_i_C(1) * t255 - r_i_i_C(2) * t256) * qJD(5)) * t258;
	t234 = t255 * t278 - t312 * t256;
	t233 = t312 * t255 + t256 * t278;
	t232 = t255 * t314 + t271 * t256;
	t231 = t271 * t255 - t256 * t314;
	t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) - t316 * t265 + t267 * t273, (-t295 - t307 + t310) * t299 + t270 * t267 + t274, (-r_i_i_C(3) * t300 - t299 * t308) * t257 + (-r_i_i_C(3) * t299 - t267 * t279) * t258 + t274, t258 * t300 - t290, t231 * r_i_i_C(1) + t232 * r_i_i_C(2), 0; -t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t265 * t273 + t316 * t267, t270 * t265 + (t275 - t310) * t298 + t276, t265 * t272 + t275 * t298 + t276, t258 * t301 + t289, -t233 * r_i_i_C(1) + t234 * r_i_i_C(2), 0; 0, t269 - t294, t269, t260 * t257, (-t256 * t303 + t288) * r_i_i_C(2) - t318 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:42
	% EndTime: 2019-10-10 11:38:43
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (709->86), mult. (745->121), div. (0->0), fcn. (606->10), ass. (0->60)
	t373 = pkin(5) + r_i_i_C(1);
	t312 = cos(pkin(10)) * pkin(4) + pkin(3);
	t318 = pkin(10) + qJ(5);
	t314 = sin(t318);
	t369 = r_i_i_C(3) + qJ(6);
	t379 = t369 * t314 + t312;
	t320 = qJ(2) + qJ(3);
	t317 = cos(t320);
	t324 = sin(qJ(1));
	t362 = qJD(1) * t324;
	t350 = t317 * t362;
	t316 = sin(t320);
	t319 = qJD(2) + qJD(3);
	t326 = cos(qJ(1));
	t365 = t319 * t326;
	t353 = t316 * t365;
	t377 = t350 + t353;
	t323 = sin(qJ(2));
	t368 = pkin(2) * qJD(2);
	t355 = t323 * t368;
	t356 = qJD(6) * t314;
	t322 = -pkin(9) - qJ(4);
	t370 = r_i_i_C(2) - t322;
	t376 = (t370 * t319 + t356) * t317 + (pkin(4) * sin(pkin(10)) + pkin(8) + pkin(7)) * qJD(1) - (t312 * t319 - qJD(4)) * t316 - t355;
	t372 = pkin(2) * t323;
	t371 = r_i_i_C(2) * t317;
	t367 = t317 * t319;
	t366 = t319 * t324;
	t364 = t324 * t314;
	t315 = cos(t318);
	t363 = t326 * t315;
	t361 = qJD(1) * t326;
	t360 = qJD(4) * t317;
	t359 = qJD(5) * t315;
	t358 = qJD(5) * t324;
	t357 = qJD(5) * t326;
	t354 = t316 * t366;
	t352 = t373 * qJD(5);
	t351 = t316 * t362;
	t348 = t314 * t358;
	t347 = t314 * t357;
	t346 = t315 * t357;
	t344 = qJD(5) * t369;
	t338 = t373 * t316 * t348 + t322 * t354 + t324 * t360 + t361 * t371;
	t337 = t317 * t363 + t364;
	t336 = t314 * t361 + t315 * t358;
	t335 = -t315 * t344 - t356;
	t334 = -t373 * t315 - t379;
	t333 = t377 * t322 + t326 * t360 + t373 * (t315 * t351 + t316 * t347) + t379 * t351;
	t325 = cos(qJ(2));
	t332 = -t315 * qJD(6) + (-t325 * pkin(2) - t312 * t317 - t370 * t316 - pkin(1)) * qJD(1);
	t331 = t334 * t316 - t317 * t322;
	t330 = t335 * t316 + (-r_i_i_C(2) * t316 + t334 * t317) * t319;
	t329 = r_i_i_C(2) * t367 + t316 * qJD(4) + t331 * t319 + (-t314 * t352 + t369 * t359 + t356) * t317;
	t328 = -t325 * t368 + t330;
	t281 = t337 * qJD(1) - t315 * t354 - t317 * t348 - t346;
	t280 = -t314 * t354 - t315 * t362 + t336 * t317 - t347;
	t279 = t377 * t315 + t317 * t347 - t336;
	t278 = t314 * t353 - t317 * t346 - t348 + (t317 * t364 + t363) * qJD(1);
	t1 = [-t369 * t280 - t373 * t281 - t376 * t324 + t332 * t326, (-t371 + t372) * t362 + t328 * t326 + t333, -r_i_i_C(2) * t350 + t330 * t326 + t333, t317 * t365 - t351, t337 * qJD(6) + t373 * t278 - t369 * t279, -t278; -t369 * t278 - t373 * t279 + t332 * t324 + t376 * t326, (t331 - t372) * t361 + t328 * t324 + t338, (-t322 * t361 + t334 * t366) * t317 + ((-r_i_i_C(2) * t319 + t335) * t324 + t334 * t361) * t316 + t338, t316 * t361 + t317 * t366, -(-t324 * t317 * t315 + t326 * t314) * qJD(6) + t369 * t281 - t373 * t280, t280; 0, t329 - t355, t329, t319 * t316, (-t316 * t344 - t373 * t367) * t314 + (t369 * t367 + (qJD(6) - t352) * t316) * t315, t314 * t367 + t316 * t359;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end