% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
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
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
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
	% StartTime: 2019-10-10 10:17:11
	% EndTime: 2019-10-10 10:17:11
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (55->22), mult. (178->38), div. (0->0), fcn. (133->6), ass. (0->17)
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t158 = sin(pkin(10));
	t159 = cos(pkin(10));
	t169 = r_i_i_C(1) * t159 - r_i_i_C(2) * t158 + pkin(2);
	t174 = r_i_i_C(3) + qJ(3);
	t175 = t169 * t160 - t174 * t162;
	t176 = t175 * qJD(2) - t160 * qJD(3);
	t161 = sin(qJ(1));
	t173 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t172 = qJD(1) * t163;
	t171 = qJD(2) * t163;
	t168 = t158 * r_i_i_C(1) + t159 * r_i_i_C(2) + pkin(7);
	t167 = -t174 * t160 - t169 * t162;
	t165 = -pkin(1) + t167;
	t1 = [t176 * t161 + (-t168 * t161 + t165 * t163) * qJD(1), (t169 * t173 - t174 * t171) * t160 + (-t174 * t173 + (-t169 * qJD(2) + qJD(3)) * t163) * t162, -t160 * t173 + t162 * t171, 0, 0, 0; -t176 * t163 + (t165 * t161 + t168 * t163) * qJD(1), -t175 * t172 + (t167 * qJD(2) + qJD(3) * t162) * t161, t161 * qJD(2) * t162 + t160 * t172, 0, 0, 0; 0, -t176, qJD(2) * t160, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:11
	% EndTime: 2019-10-10 10:17:11
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (170->42), mult. (303->71), div. (0->0), fcn. (242->8), ass. (0->37)
	t203 = cos(pkin(10)) * pkin(3) + pkin(2);
	t209 = sin(qJ(2));
	t227 = t209 * qJD(3);
	t211 = cos(qJ(2));
	t236 = r_i_i_C(3) + pkin(8) + qJ(3);
	t238 = t236 * t211;
	t242 = (-t203 * t209 + t238) * qJD(2) + t227;
	t206 = pkin(10) + qJ(4);
	t204 = sin(t206);
	t205 = cos(t206);
	t217 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
	t214 = -t217 * t209 + t238;
	t212 = cos(qJ(1));
	t228 = qJD(4) * t211;
	t221 = -qJD(1) + t228;
	t240 = t212 * t221;
	t210 = sin(qJ(1));
	t220 = qJD(1) * t211 - qJD(4);
	t232 = qJD(2) * t209;
	t237 = -t210 * t232 + t220 * t212;
	t234 = qJD(1) * t210;
	t233 = qJD(1) * t212;
	t231 = qJD(2) * t211;
	t230 = qJD(2) * t212;
	t229 = qJD(4) * t209;
	t226 = pkin(3) * sin(pkin(10)) + pkin(7);
	t224 = t236 * t209;
	t219 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
	t218 = t221 * t210;
	t216 = -t203 * t211 - pkin(1) - t224;
	t215 = t209 * t230 + t220 * t210;
	t213 = qJD(3) * t211 + t219 * t229 + (-t217 * t211 - t224) * qJD(2);
	t202 = t204 * t218 - t237 * t205;
	t201 = t237 * t204 + t205 * t218;
	t200 = t204 * t240 + t215 * t205;
	t199 = t215 * t204 - t205 * t240;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t242 * t210 + (-t226 * t210 + t216 * t212) * qJD(1), t213 * t212 - t214 * t234, -t209 * t234 + t211 * t230, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t242 * t212 + (t216 * t210 + t226 * t212) * qJD(1), t213 * t210 + t214 * t233, t209 * t233 + t210 * t231, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0, 0; 0, t214 * qJD(2) - t219 * t228 + t227, t232, (t204 * t229 - t205 * t231) * r_i_i_C(2) + (-t204 * t231 - t205 * t229) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:12
	% EndTime: 2019-10-10 10:17:12
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (335->61), mult. (543->94), div. (0->0), fcn. (458->8), ass. (0->45)
	t259 = sin(qJ(2));
	t253 = cos(pkin(10)) * pkin(3) + pkin(2);
	t256 = pkin(10) + qJ(4);
	t254 = sin(t256);
	t255 = cos(t256);
	t293 = r_i_i_C(3) + qJ(5);
	t295 = -r_i_i_C(1) - pkin(4);
	t296 = t293 * t254 - t295 * t255 + t253;
	t261 = cos(qJ(2));
	t294 = r_i_i_C(2) + pkin(8) + qJ(3);
	t298 = t294 * t261;
	t302 = t296 * t259 - t298;
	t279 = t259 * qJD(3);
	t281 = qJD(5) * t254;
	t301 = (-t253 * t259 + t298) * qJD(2) + t261 * t281 + t279;
	t299 = t281 + (t295 * t254 + t293 * t255) * qJD(4);
	t260 = sin(qJ(1));
	t291 = t260 * t261;
	t262 = cos(qJ(1));
	t290 = t262 * t255;
	t289 = qJD(1) * t260;
	t288 = qJD(1) * t262;
	t287 = qJD(2) * t259;
	t286 = qJD(2) * t261;
	t285 = qJD(2) * t262;
	t284 = qJD(4) * t259;
	t283 = qJD(4) * t260;
	t282 = qJD(4) * t262;
	t280 = t255 * qJD(5);
	t278 = pkin(3) * sin(pkin(10)) + pkin(7);
	t277 = t260 * t287;
	t276 = t254 * t283;
	t275 = t259 * t285;
	t274 = t254 * t282;
	t273 = t255 * t282;
	t272 = t294 * t259;
	t269 = t260 * t254 + t261 * t290;
	t267 = -t253 * t261 - pkin(1) - t272;
	t266 = t254 * t288 + t255 * t283;
	t263 = qJD(3) * t261 - t299 * t259 + (-t261 * t296 - t272) * qJD(2);
	t242 = t269 * qJD(1) - t255 * t277 - t261 * t276 - t273;
	t241 = -t254 * t277 - t255 * t289 + t266 * t261 - t274;
	t240 = t261 * t274 + (t261 * t289 + t275) * t255 - t266;
	t239 = t254 * t275 - t261 * t273 - t276 + (t254 * t291 + t290) * qJD(1);
	t1 = [-t262 * t280 + t295 * t242 - t293 * t241 - t301 * t260 + (-t278 * t260 + t267 * t262) * qJD(1), t263 * t262 + t302 * t289, -t259 * t289 + t261 * t285, t269 * qJD(5) - t295 * t239 - t293 * t240, -t239, 0; -t260 * t280 + t295 * t240 - t293 * t239 + t301 * t262 + (t267 * t260 + t278 * t262) * qJD(1), t263 * t260 - t288 * t302, t259 * t288 + t260 * t286, -(t262 * t254 - t255 * t291) * qJD(5) + t293 * t242 + t295 * t241, t241, 0; 0, -qJD(2) * t302 + t299 * t261 + t279, t287, (-t293 * t284 + t295 * t286) * t254 + (t293 * t286 + (t295 * qJD(4) + qJD(5)) * t259) * t255, t254 * t286 + t255 * t284, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:12
	% EndTime: 2019-10-10 10:17:13
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (753->95), mult. (1191->147), div. (0->0), fcn. (1113->10), ass. (0->63)
	t302 = sin(qJ(2));
	t295 = cos(pkin(10)) * pkin(3) + pkin(2);
	t298 = pkin(10) + qJ(4);
	t296 = sin(t298);
	t297 = cos(t298);
	t304 = cos(qJ(6));
	t301 = sin(qJ(6));
	t346 = -pkin(4) - pkin(5);
	t323 = t301 * r_i_i_C(2) + t346;
	t312 = t304 * r_i_i_C(1) - t323;
	t324 = -t301 * r_i_i_C(1) - qJ(5);
	t313 = t304 * r_i_i_C(2) - t324;
	t347 = t313 * t296 + t312 * t297 + t295;
	t305 = cos(qJ(2));
	t331 = -r_i_i_C(3) - pkin(9) + pkin(8) + qJ(3);
	t353 = t331 * t305;
	t358 = t347 * t302 - t353;
	t357 = (-t302 * t295 + t353) * qJD(2) + t302 * qJD(3);
	t306 = cos(qJ(1));
	t342 = t306 * t297;
	t303 = sin(qJ(1));
	t344 = t303 * t305;
	t279 = t296 * t344 + t342;
	t343 = t306 * t296;
	t280 = t297 * t344 - t343;
	t318 = t279 * t301 + t280 * t304;
	t319 = t279 * t304 - t280 * t301;
	t355 = (t318 * r_i_i_C(1) + t319 * r_i_i_C(2)) * qJD(6);
	t314 = t296 * t301 + t297 * t304;
	t315 = t296 * t304 - t297 * t301;
	t354 = (t312 * t296 - t313 * t297) * qJD(4) - (t315 * r_i_i_C(1) - t314 * r_i_i_C(2)) * qJD(6) - t331 * qJD(2) - t296 * qJD(5);
	t351 = (-qJD(4) + qJD(6)) * t302;
	t335 = qJD(4) * t303;
	t340 = qJD(1) * t306;
	t311 = t296 * t340 + t297 * t335;
	t334 = qJD(4) * t306;
	t326 = t296 * t334;
	t339 = qJD(2) * t302;
	t329 = t303 * t339;
	t341 = qJD(1) * t303;
	t277 = -t296 * t329 - t297 * t341 + t311 * t305 - t326;
	t274 = t277 * t304;
	t338 = qJD(2) * t305;
	t337 = qJD(2) * t306;
	t336 = qJD(4) * t302;
	t330 = pkin(3) * sin(pkin(10)) + pkin(7);
	t328 = t296 * t335;
	t327 = t302 * t337;
	t325 = t297 * t334;
	t320 = -(t314 * t351 - t315 * t338) * r_i_i_C(1) - (t314 * t338 + t315 * t351) * r_i_i_C(2);
	t281 = -t303 * t297 + t305 * t343;
	t282 = t303 * t296 + t305 * t342;
	t317 = t281 * t304 - t282 * t301;
	t316 = t281 * t301 + t282 * t304;
	t310 = -t295 * t305 - t331 * t302 - pkin(1);
	t308 = -qJD(2) * t347 + qJD(3);
	t307 = t354 * t302 + t308 * t305;
	t278 = t282 * qJD(1) - t297 * t329 - t305 * t328 - t325;
	t276 = t305 * t326 + (t305 * t341 + t327) * t297 - t311;
	t275 = t279 * qJD(1) + t296 * t327 - t305 * t325 - t328;
	t271 = t317 * qJD(6) - t275 * t301 - t276 * t304;
	t270 = -t316 * qJD(6) - t275 * t304 + t276 * t301;
	t1 = [-t274 * r_i_i_C(2) - t279 * qJD(5) + t324 * t277 - t312 * t278 + (-t319 * r_i_i_C(1) + t318 * r_i_i_C(2)) * qJD(6) - t357 * t303 + (-t330 * t303 + t310 * t306) * qJD(1), t307 * t306 + t358 * t341, -t302 * t341 + t305 * t337, t282 * qJD(5) - t313 * t276 + t312 * t275 + (t316 * r_i_i_C(1) + t317 * r_i_i_C(2)) * qJD(6), -t275, t270 * r_i_i_C(1) - t271 * r_i_i_C(2); t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t275 * qJ(5) + t281 * qJD(5) + t346 * t276 + t357 * t306 + (t310 * t303 + t330 * t306) * qJD(1), t307 * t303 - t358 * t340, t302 * t340 + t303 * t338, -t274 * r_i_i_C(1) + t280 * qJD(5) + t323 * t277 + t313 * t278 + t355, t277, (-t278 * t301 + t274) * r_i_i_C(1) + (-t277 * t301 - t278 * t304) * r_i_i_C(2) - t355; 0, t308 * t302 - t354 * t305, t339, (-qJ(5) * t336 + t346 * t338) * t296 + (qJ(5) * t338 + (t346 * qJD(4) + qJD(5)) * t302) * t297 - t320, t296 * t338 + t297 * t336, t320;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end