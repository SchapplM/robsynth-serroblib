% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
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
	% StartTime: 2019-10-10 12:40:25
	% EndTime: 2019-10-10 12:40:25
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
	% StartTime: 2019-10-10 12:40:26
	% EndTime: 2019-10-10 12:40:26
	% DurationCPUTime: 0.21s
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
	% StartTime: 2019-10-10 12:40:26
	% EndTime: 2019-10-10 12:40:26
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (261->56), mult. (420->90), div. (0->0), fcn. (332->8), ass. (0->52)
	t235 = cos(qJ(3));
	t227 = t235 * pkin(3) + pkin(2);
	t233 = sin(qJ(2));
	t236 = cos(qJ(2));
	t267 = r_i_i_C(3) + pkin(9) + pkin(8);
	t251 = t267 * t236;
	t232 = sin(qJ(3));
	t266 = pkin(3) * qJD(3);
	t256 = t232 * t266;
	t276 = (-t227 * t233 + t251) * qJD(2) - t236 * t256;
	t237 = cos(qJ(1));
	t230 = qJD(3) + qJD(4);
	t250 = t230 * t236 - qJD(1);
	t274 = t237 * t250;
	t231 = qJ(3) + qJ(4);
	t228 = sin(t231);
	t229 = cos(t231);
	t268 = r_i_i_C(2) * t229;
	t246 = r_i_i_C(1) * t228 + t268;
	t273 = -t246 * t230 - t256;
	t260 = qJD(1) * t236;
	t249 = -t230 + t260;
	t234 = sin(qJ(1));
	t258 = qJD(2) * t233;
	t253 = t234 * t258;
	t272 = t249 * t237 - t253;
	t271 = pkin(3) * t232;
	t270 = r_i_i_C(1) * t229;
	t269 = r_i_i_C(2) * t228;
	t264 = t230 * t233;
	t252 = t237 * t258;
	t240 = t249 * t234 + t252;
	t222 = t240 * t228 - t229 * t274;
	t223 = t228 * t274 + t240 * t229;
	t263 = t222 * r_i_i_C(1) + t223 * r_i_i_C(2);
	t245 = t250 * t234;
	t224 = t272 * t228 + t229 * t245;
	t225 = t228 * t245 - t272 * t229;
	t262 = -t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
	t261 = qJD(1) * t234;
	t259 = qJD(1) * t237;
	t257 = qJD(2) * t236;
	t255 = t235 * t266;
	t254 = pkin(7) + t271;
	t247 = -qJD(3) + t260;
	t244 = (-qJD(3) * t236 + qJD(1)) * t235;
	t243 = t227 - t269 + t270;
	t242 = -t227 * t236 - t267 * t233 - pkin(1);
	t241 = qJD(2) * t243;
	t239 = -t267 * qJD(2) - t273;
	t226 = t264 * t269;
	t1 = [t237 * t255 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) - t276 * t234 + (-t254 * t234 + t242 * t237) * qJD(1), (-t237 * t241 - t267 * t261) * t236 + (t239 * t237 + t243 * t261) * t233, (t237 * t244 + (t247 * t234 + t252) * t232) * pkin(3) + t263, t263, 0, 0; t234 * t255 - t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t276 * t237 + (t242 * t234 + t254 * t237) * qJD(1), (-t234 * t241 + t267 * t259) * t236 + (t239 * t234 - t243 * t259) * t233, (t234 * t244 + (-t247 * t237 + t253) * t232) * pkin(3) + t262, t262, 0, 0; 0, t273 * t236 + (-t243 * t233 + t251) * qJD(2), t226 + (-t230 * t270 - t255) * t233 + (-t246 - t271) * t257, -t257 * t268 + t226 + (-t228 * t257 - t229 * t264) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:26
	% EndTime: 2019-10-10 12:40:27
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (462->65), mult. (522->93), div. (0->0), fcn. (410->10), ass. (0->58)
	t248 = qJ(3) + qJ(4);
	t243 = sin(t248);
	t249 = sin(qJ(3));
	t284 = pkin(3) * qJD(3);
	t247 = qJD(3) + qJD(4);
	t289 = pkin(4) * t247;
	t233 = -t243 * t289 - t249 * t284;
	t244 = cos(t248);
	t241 = pkin(4) * t244;
	t252 = cos(qJ(3));
	t238 = t252 * pkin(3) + t241;
	t236 = pkin(2) + t238;
	t290 = pkin(4) * t243;
	t237 = t249 * pkin(3) + t290;
	t250 = sin(qJ(2));
	t253 = cos(qJ(2));
	t273 = t250 * qJD(5);
	t285 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
	t294 = t285 * t253;
	t297 = (-t236 * t250 + t294) * qJD(2) + (pkin(7) + t237) * qJD(1) + t253 * t233 + t273;
	t242 = pkin(11) + t248;
	t239 = sin(t242);
	t287 = r_i_i_C(2) * t239;
	t240 = cos(t242);
	t288 = r_i_i_C(1) * t240;
	t261 = t236 - t287 + t288;
	t256 = -t261 * t250 + t294;
	t251 = sin(qJ(1));
	t268 = t247 * t253 - qJD(1);
	t263 = t268 * t251;
	t254 = cos(qJ(1));
	t264 = t268 * t254;
	t266 = r_i_i_C(1) * t239 + r_i_i_C(2) * t240;
	t293 = t266 * t247 - t233;
	t278 = qJD(1) * t253;
	t267 = -t247 + t278;
	t276 = qJD(2) * t250;
	t292 = -t251 * t276 + t267 * t254;
	t282 = t247 * t250;
	t274 = qJD(2) * t254;
	t258 = t250 * t274 + t267 * t251;
	t229 = t258 * t239 - t240 * t264;
	t230 = t239 * t264 + t258 * t240;
	t281 = t229 * r_i_i_C(1) + t230 * r_i_i_C(2);
	t231 = t292 * t239 + t240 * t263;
	t232 = t239 * t263 - t240 * t292;
	t280 = -t231 * r_i_i_C(1) + t232 * r_i_i_C(2);
	t279 = qJD(1) * t251;
	t277 = qJD(1) * t254;
	t275 = qJD(2) * t253;
	t271 = t285 * t250;
	t265 = t237 * t278 + t233;
	t234 = t244 * t289 + t252 * t284;
	t260 = qJD(1) * t238 - t234 * t253 + t237 * t276;
	t257 = t234 + (-t236 * t253 - pkin(1) - t271) * qJD(1);
	t255 = qJD(5) * t253 + t293 * t250 + (-t261 * t253 - t271) * qJD(2);
	t235 = t282 * t287;
	t1 = [t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t297 * t251 + t257 * t254, t255 * t254 - t256 * t279, t265 * t251 + t260 * t254 + t281, (t258 * t243 - t244 * t264) * pkin(4) + t281, -t250 * t279 + t253 * t274, 0; -t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t257 * t251 + t297 * t254, t255 * t251 + t256 * t277, t260 * t251 - t265 * t254 + t280, (-t243 * t292 - t244 * t263) * pkin(4) + t280, t250 * t277 + t251 * t275, 0; 0, t256 * qJD(2) - t293 * t253 + t273, t235 + (-t247 * t288 - t234) * t250 + (-t237 - t266) * t275, t235 + (-t288 - t241) * t282 + (-t266 - t290) * t275, t276, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:40:26
	% EndTime: 2019-10-10 12:40:27
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (845->75), mult. (649->100), div. (0->0), fcn. (512->12), ass. (0->63)
	t260 = qJ(3) + qJ(4);
	t254 = pkin(11) + t260;
	t244 = -pkin(4) * sin(t260) - pkin(5) * sin(t254);
	t259 = qJD(3) + qJD(4);
	t238 = t244 * t259;
	t261 = sin(qJ(3));
	t298 = pkin(3) * qJD(3);
	t236 = -t261 * t298 + t238;
	t264 = cos(qJ(3));
	t307 = pkin(5) * cos(t254) + pkin(4) * cos(t260);
	t242 = t264 * pkin(3) + t307;
	t240 = pkin(2) + t242;
	t241 = t261 * pkin(3) - t244;
	t262 = sin(qJ(2));
	t265 = cos(qJ(2));
	t287 = t262 * qJD(5);
	t299 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9) + pkin(8);
	t308 = t299 * t265;
	t312 = (-t240 * t262 + t308) * qJD(2) + (pkin(7) + t241) * qJD(1) + t265 * t236 + t287;
	t252 = qJ(6) + t254;
	t247 = sin(t252);
	t302 = r_i_i_C(2) * t247;
	t248 = cos(t252);
	t303 = r_i_i_C(1) * t248;
	t273 = t240 - t302 + t303;
	t268 = -t273 * t262 + t308;
	t266 = cos(qJ(1));
	t253 = qJD(6) + t259;
	t281 = t253 * t265 - qJD(1);
	t310 = t266 * t281;
	t239 = t307 * t259;
	t301 = r_i_i_C(2) * t248;
	t278 = r_i_i_C(1) * t247 + t301;
	t306 = t278 * t253 - t236;
	t263 = sin(qJ(1));
	t292 = qJD(1) * t265;
	t280 = -t253 + t292;
	t290 = qJD(2) * t262;
	t305 = -t263 * t290 + t280 * t266;
	t296 = t253 * t262;
	t288 = qJD(2) * t266;
	t270 = t262 * t288 + t280 * t263;
	t232 = t270 * t247 - t248 * t310;
	t233 = t247 * t310 + t270 * t248;
	t295 = t232 * r_i_i_C(1) + t233 * r_i_i_C(2);
	t275 = t281 * t263;
	t234 = t305 * t247 + t248 * t275;
	t235 = t247 * t275 - t305 * t248;
	t294 = -t234 * r_i_i_C(1) + t235 * r_i_i_C(2);
	t293 = qJD(1) * t263;
	t291 = qJD(1) * t266;
	t289 = qJD(2) * t265;
	t286 = t253 * t303;
	t284 = t299 * t262;
	t277 = t241 * t292 + t236;
	t276 = t244 * t292 - t238;
	t237 = t264 * t298 + t239;
	t272 = qJD(1) * t242 - t237 * t265 + t241 * t290;
	t271 = qJD(1) * t307 - t239 * t265 - t244 * t290;
	t269 = t237 + (-t240 * t265 - pkin(1) - t284) * qJD(1);
	t267 = qJD(5) * t265 + t306 * t262 + (-t273 * t265 - t284) * qJD(2);
	t243 = t296 * t302;
	t1 = [t235 * r_i_i_C(1) + t234 * r_i_i_C(2) - t312 * t263 + t269 * t266, t267 * t266 - t268 * t293, t277 * t263 + t272 * t266 + t295, -t276 * t263 + t271 * t266 + t295, -t262 * t293 + t265 * t288, t295; -t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t269 * t263 + t312 * t266, t267 * t263 + t268 * t291, t272 * t263 - t277 * t266 + t294, t271 * t263 + t276 * t266 + t294, t262 * t291 + t263 * t289, t294; 0, t268 * qJD(2) - t306 * t265 + t287, t243 + (-t237 - t286) * t262 + (-t241 - t278) * t289, t243 + (-t239 - t286) * t262 + (t244 - t278) * t289, t290, -t289 * t301 + t243 + (-t247 * t289 - t248 * t296) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end