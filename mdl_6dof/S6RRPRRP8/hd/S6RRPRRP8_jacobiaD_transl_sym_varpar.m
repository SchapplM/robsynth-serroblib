% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP8
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
% Datum: 2019-10-10 10:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
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
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
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
	% StartTime: 2019-10-10 10:41:04
	% EndTime: 2019-10-10 10:41:05
	% DurationCPUTime: 0.21s
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
	% StartTime: 2019-10-10 10:41:05
	% EndTime: 2019-10-10 10:41:05
	% DurationCPUTime: 0.23s
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
	t226 = sin(pkin(10)) * pkin(3) + pkin(7);
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
	% StartTime: 2019-10-10 10:41:05
	% EndTime: 2019-10-10 10:41:05
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (402->59), mult. (449->88), div. (0->0), fcn. (361->10), ass. (0->55)
	t238 = pkin(10) + qJ(4);
	t235 = cos(t238);
	t230 = pkin(4) * t235 + cos(pkin(10)) * pkin(3) + pkin(2);
	t240 = sin(qJ(2));
	t242 = cos(qJ(2));
	t234 = sin(t238);
	t274 = pkin(4) * qJD(4);
	t262 = t234 * t274;
	t263 = t240 * qJD(3);
	t280 = pkin(4) * t234;
	t275 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
	t284 = t275 * t242;
	t288 = (-t230 * t240 + t284) * qJD(2) + (pkin(7) + t280 + sin(pkin(10)) * pkin(3)) * qJD(1) - t242 * t262 + t263;
	t236 = qJ(5) + t238;
	t232 = sin(t236);
	t278 = r_i_i_C(2) * t232;
	t233 = cos(t236);
	t279 = r_i_i_C(1) * t233;
	t249 = t230 - t278 + t279;
	t246 = -t249 * t240 + t284;
	t243 = cos(qJ(1));
	t239 = qJD(4) + qJD(5);
	t255 = t239 * t242 - qJD(1);
	t286 = t243 * t255;
	t277 = r_i_i_C(2) * t233;
	t252 = r_i_i_C(1) * t232 + t277;
	t283 = t252 * t239 + t262;
	t268 = qJD(1) * t242;
	t254 = -t239 + t268;
	t241 = sin(qJ(1));
	t266 = qJD(2) * t240;
	t260 = t241 * t266;
	t282 = t254 * t243 - t260;
	t272 = t239 * t240;
	t264 = qJD(2) * t243;
	t259 = t240 * t264;
	t247 = t254 * t241 + t259;
	t225 = t247 * t232 - t233 * t286;
	t226 = t232 * t286 + t247 * t233;
	t271 = t225 * r_i_i_C(1) + t226 * r_i_i_C(2);
	t251 = t255 * t241;
	t227 = t282 * t232 + t233 * t251;
	t228 = t232 * t251 - t282 * t233;
	t270 = -t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
	t269 = qJD(1) * t241;
	t267 = qJD(1) * t243;
	t265 = qJD(2) * t242;
	t261 = t235 * t274;
	t258 = t275 * t240;
	t253 = -qJD(4) + t268;
	t250 = t235 * (-qJD(4) * t242 + qJD(1));
	t245 = t261 + (-t230 * t242 - pkin(1) - t258) * qJD(1);
	t244 = qJD(3) * t242 + t283 * t240 + (-t249 * t242 - t258) * qJD(2);
	t229 = t272 * t278;
	t1 = [t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t288 * t241 + t245 * t243, t244 * t243 - t246 * t269, -t240 * t269 + t242 * t264, (t243 * t250 + (t253 * t241 + t259) * t234) * pkin(4) + t271, t271, 0; -t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t245 * t241 + t288 * t243, t244 * t241 + t246 * t267, t240 * t267 + t241 * t265, (t241 * t250 + (-t253 * t243 + t260) * t234) * pkin(4) + t270, t270, 0; 0, t246 * qJD(2) - t283 * t242 + t263, t266, t229 + (-t239 * t279 - t261) * t240 + (-t252 - t280) * t265, -t265 * t277 + t229 + (-t232 * t265 - t233 * t272) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:05
	% EndTime: 2019-10-10 10:41:06
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (777->73), mult. (767->104), div. (0->0), fcn. (649->10), ass. (0->57)
	t304 = pkin(10) + qJ(4);
	t302 = qJ(5) + t304;
	t299 = cos(t302);
	t348 = r_i_i_C(3) + qJ(6);
	t362 = t348 * t299;
	t305 = qJD(4) + qJD(5);
	t298 = sin(t302);
	t300 = sin(t304);
	t347 = pkin(4) * qJD(4);
	t349 = r_i_i_C(2) + pkin(9) + pkin(8) + qJ(3);
	t314 = t349 * qJD(2) + qJD(6) * t298 - t300 * t347;
	t352 = -r_i_i_C(1) - pkin(5);
	t333 = t352 * t298;
	t361 = (t333 + t362) * t305 + t314;
	t301 = cos(t304);
	t291 = pkin(4) * t301 + cos(pkin(10)) * pkin(3) + pkin(2);
	t306 = sin(qJ(2));
	t308 = cos(qJ(2));
	t351 = pkin(4) * t300;
	t360 = t314 * t308 + (pkin(7) + t351 + sin(pkin(10)) * pkin(3)) * qJD(1) - (qJD(2) * t291 - qJD(3)) * t306;
	t318 = -t348 * t298 + t352 * t299;
	t315 = -t291 + t318;
	t358 = t315 * t306 + t349 * t308;
	t307 = sin(qJ(1));
	t309 = cos(qJ(1));
	t342 = t309 * t299;
	t320 = t307 * t298 + t308 * t342;
	t344 = t307 * t308;
	t354 = t298 * t344 + t342;
	t346 = t305 * t306;
	t343 = t309 * t298;
	t341 = qJD(1) * t307;
	t340 = qJD(1) * t308;
	t339 = qJD(1) * t309;
	t338 = qJD(2) * t306;
	t337 = qJD(2) * t308;
	t336 = qJD(2) * t309;
	t335 = t299 * qJD(6);
	t334 = t301 * t347;
	t332 = t299 * t346;
	t330 = t305 * t343;
	t328 = t306 * t335 + t337 * t362;
	t326 = t307 * t338;
	t325 = t306 * t336;
	t324 = -qJD(4) + t340;
	t322 = t301 * (-qJD(4) * t308 + qJD(1));
	t319 = t307 * t305 * t299 + t298 * t339;
	t277 = t354 * qJD(1) + t298 * t325 - t320 * t305;
	t278 = t308 * t330 + (t307 * t340 + t325) * t299 - t319;
	t317 = qJD(6) * t320 - t352 * t277 - t348 * t278;
	t279 = -t298 * t326 - t299 * t341 + t319 * t308 - t330;
	t280 = t320 * qJD(1) - t299 * t326 - t354 * t305;
	t316 = -(-t299 * t344 + t343) * qJD(6) + t348 * t280 + t352 * t279;
	t312 = t315 * qJD(2) + qJD(3);
	t311 = t334 - t335 + (-t291 * t308 - t349 * t306 - pkin(1)) * qJD(1);
	t310 = -t361 * t306 + t312 * t308;
	t1 = [-t348 * t279 + t352 * t280 - t360 * t307 + t311 * t309, t310 * t309 - t358 * t341, -t306 * t341 + t308 * t336, (t309 * t322 + (t324 * t307 + t325) * t300) * pkin(4) + t317, t317, -t277; -t348 * t277 + t352 * t278 + t311 * t307 + t360 * t309, t310 * t307 + t358 * t339, t306 * t339 + t307 * t337, (t307 * t322 + (-t324 * t309 + t326) * t300) * pkin(4) + t316, t316, t279; 0, t312 * t306 + t361 * t308, t338, (t333 - t351) * t337 + (t318 * t305 - t334) * t306 + t328, t352 * t332 + (t352 * t337 - t348 * t346) * t298 + t328, t298 * t337 + t332;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end