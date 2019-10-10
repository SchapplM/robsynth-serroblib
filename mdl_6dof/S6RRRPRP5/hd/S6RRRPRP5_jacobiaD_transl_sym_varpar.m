% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP5
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
% Datum: 2019-10-10 11:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:11
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:11
	% EndTime: 2019-10-10 11:42:11
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
	% StartTime: 2019-10-10 11:42:11
	% EndTime: 2019-10-10 11:42:11
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 11:42:12
	% EndTime: 2019-10-10 11:42:12
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t223 = pkin(8) + r_i_i_C(3);
	t213 = t223 * t202;
	t224 = -pkin(2) * t199 + t213;
	t201 = cos(qJ(3));
	t203 = cos(qJ(1));
	t221 = t201 * t203;
	t200 = sin(qJ(1));
	t220 = qJD(1) * t200;
	t219 = qJD(1) * t203;
	t218 = qJD(2) * t200;
	t217 = qJD(2) * t202;
	t216 = qJD(2) * t203;
	t215 = qJD(3) * t199;
	t214 = qJD(3) * t202;
	t212 = -qJD(1) + t214;
	t211 = qJD(1) * t202 - qJD(3);
	t198 = sin(qJ(3));
	t210 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201;
	t209 = r_i_i_C(1) * t201 - r_i_i_C(2) * t198 + pkin(2);
	t208 = t212 * t198;
	t207 = -pkin(2) * t202 - t223 * t199 - pkin(1);
	t206 = qJD(2) * t209;
	t205 = t199 * t216 + t211 * t200;
	t204 = -t223 * qJD(2) + t210 * qJD(3);
	t197 = -t211 * t221 + (qJD(2) * t199 * t201 + t208) * t200;
	t196 = t212 * t201 * t200 + (-t199 * t218 + t211 * t203) * t198;
	t195 = t205 * t201 + t203 * t208;
	t194 = t205 * t198 - t212 * t221;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(7) * t200 + t207 * t203) * qJD(1), (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0, 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(7) * t203 + t207 * t200) * qJD(1), (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0, 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2), (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:12
	% EndTime: 2019-10-10 11:42:12
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (189->46), mult. (370->74), div. (0->0), fcn. (290->8), ass. (0->40)
	t213 = cos(qJ(3));
	t242 = t213 * pkin(3);
	t205 = pkin(2) + t242;
	t211 = sin(qJ(2));
	t214 = cos(qJ(2));
	t229 = qJD(3) * t214 - qJD(1);
	t234 = t211 * qJD(4);
	t210 = sin(qJ(3));
	t243 = pkin(3) * t210;
	t241 = r_i_i_C(3) + qJ(4) + pkin(8);
	t246 = t241 * t214;
	t249 = pkin(7) * qJD(1) + (-t205 * t211 + t246) * qJD(2) - t229 * t243 + t234;
	t208 = qJ(3) + pkin(10);
	t206 = sin(t208);
	t207 = cos(t208);
	t227 = r_i_i_C(1) * t207 - r_i_i_C(2) * t206;
	t223 = t205 + t227;
	t218 = -t223 * t211 + t246;
	t212 = sin(qJ(1));
	t224 = t229 * t212;
	t215 = cos(qJ(1));
	t225 = t229 * t215;
	t228 = qJD(1) * t214 - qJD(3);
	t237 = qJD(2) * t212;
	t245 = -t211 * t237 + t228 * t215;
	t239 = qJD(1) * t212;
	t238 = qJD(1) * t215;
	t236 = qJD(2) * t215;
	t235 = qJD(3) * t211;
	t232 = t241 * t211;
	t222 = r_i_i_C(1) * t206 + r_i_i_C(2) * t207 + t243;
	t221 = t222 * t214;
	t219 = t211 * t236 + t228 * t212;
	t217 = qJD(3) * t242 + (-t205 * t214 - pkin(1) - t232) * qJD(1);
	t216 = qJD(4) * t214 + t222 * t235 + (-t223 * t214 - t232) * qJD(2);
	t204 = t206 * t224 - t207 * t245;
	t203 = t245 * t206 + t207 * t224;
	t202 = t206 * t225 + t219 * t207;
	t201 = t219 * t206 - t207 * t225;
	t1 = [t204 * r_i_i_C(1) + t203 * r_i_i_C(2) - t249 * t212 + t217 * t215, t216 * t215 - t218 * t239, t201 * r_i_i_C(1) + t202 * r_i_i_C(2) + (t219 * t210 - t213 * t225) * pkin(3), -t211 * t239 + t214 * t236, 0, 0; -t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t217 * t212 + t249 * t215, t216 * t212 + t218 * t238, -t203 * r_i_i_C(1) + t204 * r_i_i_C(2) + (-t210 * t245 - t213 * t224) * pkin(3), t211 * t238 + t214 * t237, 0, 0; 0, t218 * qJD(2) - qJD(3) * t221 + t234, (-t227 - t242) * t235 - qJD(2) * t221, qJD(2) * t211, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:12
	% EndTime: 2019-10-10 11:42:12
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (421->58), mult. (480->85), div. (0->0), fcn. (380->10), ass. (0->51)
	t244 = qJ(3) + pkin(10);
	t233 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t244);
	t229 = t233 * qJD(3);
	t234 = pkin(4) * cos(t244) + cos(qJ(3)) * pkin(3);
	t232 = pkin(2) + t234;
	t246 = sin(qJ(2));
	t249 = cos(qJ(2));
	t267 = t246 * qJD(4);
	t278 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
	t286 = t278 * t249;
	t290 = (-t232 * t246 + t286) * qJD(2) + (pkin(7) + t233) * qJD(1) - t249 * t229 + t267;
	t240 = qJ(5) + t244;
	t236 = sin(t240);
	t281 = r_i_i_C(2) * t236;
	t237 = cos(t240);
	t282 = r_i_i_C(1) * t237;
	t256 = t232 - t281 + t282;
	t252 = -t256 * t246 + t286;
	t250 = cos(qJ(1));
	t243 = qJD(3) + qJD(5);
	t262 = t243 * t249 - qJD(1);
	t288 = t250 * t262;
	t280 = r_i_i_C(2) * t237;
	t260 = r_i_i_C(1) * t236 + t280;
	t285 = t260 * t243 + t229;
	t247 = sin(qJ(1));
	t272 = qJD(1) * t249;
	t261 = -t243 + t272;
	t270 = qJD(2) * t246;
	t284 = -t247 * t270 + t261 * t250;
	t276 = t243 * t246;
	t268 = qJD(2) * t250;
	t254 = t246 * t268 + t261 * t247;
	t225 = t254 * t236 - t237 * t288;
	t226 = t236 * t288 + t254 * t237;
	t275 = t225 * r_i_i_C(1) + t226 * r_i_i_C(2);
	t258 = t262 * t247;
	t227 = t284 * t236 + t237 * t258;
	t228 = t236 * t258 - t284 * t237;
	t274 = -t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
	t273 = qJD(1) * t247;
	t271 = qJD(1) * t250;
	t269 = qJD(2) * t249;
	t265 = t278 * t246;
	t259 = t233 * t272 - t229;
	t230 = t234 * qJD(3);
	t255 = qJD(1) * t234 - t230 * t249 + t233 * t270;
	t253 = t230 + (-t232 * t249 - pkin(1) - t265) * qJD(1);
	t251 = qJD(4) * t249 + t285 * t246 + (-t256 * t249 - t265) * qJD(2);
	t231 = t276 * t281;
	t1 = [t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t290 * t247 + t253 * t250, t251 * t250 - t252 * t273, t259 * t247 + t255 * t250 + t275, -t246 * t273 + t249 * t268, t275, 0; -t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t253 * t247 + t290 * t250, t251 * t247 + t252 * t271, t255 * t247 - t259 * t250 + t274, t246 * t271 + t247 * t269, t274, 0; 0, t252 * qJD(2) - t285 * t249 + t267, t231 + (-t243 * t282 - t230) * t246 + (-t233 - t260) * t269, t270, -t269 * t280 + t231 + (-t236 * t269 - t237 * t276) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:42:13
	% EndTime: 2019-10-10 11:42:13
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (796->72), mult. (798->101), div. (0->0), fcn. (668->10), ass. (0->56)
	t311 = qJ(3) + pkin(10);
	t307 = qJ(5) + t311;
	t304 = cos(t307);
	t354 = r_i_i_C(3) + qJ(6);
	t367 = t354 * t304;
	t310 = qJD(3) + qJD(5);
	t297 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t311);
	t286 = t297 * qJD(3);
	t303 = sin(t307);
	t355 = r_i_i_C(2) + pkin(9) + qJ(4) + pkin(8);
	t325 = t355 * qJD(2) + qJD(6) * t303 - t286;
	t357 = pkin(5) + r_i_i_C(1);
	t341 = t357 * t303;
	t366 = (t341 - t367) * t310 - t325;
	t298 = pkin(4) * cos(t311) + cos(qJ(3)) * pkin(3);
	t294 = pkin(2) + t298;
	t313 = sin(qJ(2));
	t316 = cos(qJ(2));
	t365 = t325 * t316 + (pkin(7) + t297) * qJD(1) - (qJD(2) * t294 - qJD(4)) * t313;
	t327 = -t354 * t303 - t357 * t304;
	t322 = -t294 + t327;
	t363 = t322 * t313 + t355 * t316;
	t314 = sin(qJ(1));
	t317 = cos(qJ(1));
	t349 = t317 * t304;
	t329 = t314 * t303 + t316 * t349;
	t351 = t314 * t316;
	t359 = t303 * t351 + t349;
	t353 = t310 * t313;
	t350 = t317 * t303;
	t348 = qJD(1) * t314;
	t347 = qJD(1) * t316;
	t346 = qJD(1) * t317;
	t345 = qJD(2) * t313;
	t344 = qJD(2) * t316;
	t343 = qJD(2) * t317;
	t342 = t304 * qJD(6);
	t340 = t304 * t353;
	t338 = t310 * t350;
	t336 = t313 * t342 + t344 * t367;
	t334 = t314 * t345;
	t333 = t313 * t343;
	t330 = t297 * t347 - t286;
	t328 = t314 * t310 * t304 + t303 * t346;
	t278 = t359 * qJD(1) + t303 * t333 - t329 * t310;
	t279 = t316 * t338 + (t314 * t347 + t333) * t304 - t328;
	t324 = qJD(6) * t329 + t357 * t278 - t354 * t279;
	t280 = -t303 * t334 - t304 * t348 + t328 * t316 - t338;
	t281 = t329 * qJD(1) - t304 * t334 - t359 * t310;
	t323 = -(-t304 * t351 + t350) * qJD(6) + t354 * t281 - t357 * t280;
	t287 = t298 * qJD(3);
	t321 = qJD(1) * t298 - t287 * t316 + t297 * t345;
	t320 = t322 * qJD(2) + qJD(4);
	t319 = -t342 + t287 + (-t294 * t316 - t355 * t313 - pkin(1)) * qJD(1);
	t318 = t366 * t313 + t320 * t316;
	t1 = [-t354 * t280 - t357 * t281 - t365 * t314 + t319 * t317, t318 * t317 - t363 * t348, t330 * t314 + t321 * t317 + t324, -t313 * t348 + t316 * t343, t324, -t278; -t354 * t278 - t357 * t279 + t319 * t314 + t365 * t317, t318 * t314 + t363 * t346, t321 * t314 - t330 * t317 + t323, t313 * t346 + t314 * t344, t323, t280; 0, t320 * t313 - t366 * t316, (-t297 - t341) * t344 + (t327 * t310 - t287) * t313 + t336, t345, -t357 * t340 + (-t357 * t344 - t354 * t353) * t303 + t336, t303 * t344 + t340;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end