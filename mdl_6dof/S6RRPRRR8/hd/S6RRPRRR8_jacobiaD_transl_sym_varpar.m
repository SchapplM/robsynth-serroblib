% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:01
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:01:24
	% EndTime: 2019-10-10 11:01:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:01:25
	% EndTime: 2019-10-10 11:01:25
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
	% StartTime: 2019-10-10 11:01:25
	% EndTime: 2019-10-10 11:01:25
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-10 11:01:26
	% EndTime: 2019-10-10 11:01:26
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (55->22), mult. (178->38), div. (0->0), fcn. (133->6), ass. (0->17)
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t158 = sin(pkin(11));
	t159 = cos(pkin(11));
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
	% StartTime: 2019-10-10 11:01:26
	% EndTime: 2019-10-10 11:01:26
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (170->42), mult. (303->71), div. (0->0), fcn. (242->8), ass. (0->37)
	t203 = cos(pkin(11)) * pkin(3) + pkin(2);
	t209 = sin(qJ(2));
	t227 = t209 * qJD(3);
	t211 = cos(qJ(2));
	t236 = r_i_i_C(3) + pkin(8) + qJ(3);
	t238 = t236 * t211;
	t242 = (-t203 * t209 + t238) * qJD(2) + t227;
	t206 = pkin(11) + qJ(4);
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
	t226 = sin(pkin(11)) * pkin(3) + pkin(7);
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
	% StartTime: 2019-10-10 11:01:26
	% EndTime: 2019-10-10 11:01:26
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (402->59), mult. (449->88), div. (0->0), fcn. (361->10), ass. (0->55)
	t238 = pkin(11) + qJ(4);
	t235 = cos(t238);
	t230 = pkin(4) * t235 + cos(pkin(11)) * pkin(3) + pkin(2);
	t240 = sin(qJ(2));
	t242 = cos(qJ(2));
	t234 = sin(t238);
	t274 = pkin(4) * qJD(4);
	t262 = t234 * t274;
	t263 = t240 * qJD(3);
	t280 = pkin(4) * t234;
	t275 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
	t284 = t275 * t242;
	t288 = (-t230 * t240 + t284) * qJD(2) + (pkin(7) + t280 + sin(pkin(11)) * pkin(3)) * qJD(1) - t242 * t262 + t263;
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
	% StartTime: 2019-10-10 11:01:26
	% EndTime: 2019-10-10 11:01:26
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (785->76), mult. (601->103), div. (0->0), fcn. (481->12), ass. (0->66)
	t256 = pkin(11) + qJ(4);
	t254 = qJ(5) + t256;
	t249 = cos(t254);
	t252 = cos(t256);
	t243 = pkin(4) * t252 + pkin(5) * t249;
	t237 = cos(pkin(11)) * pkin(3) + pkin(2) + t243;
	t248 = sin(t254);
	t251 = sin(t256);
	t294 = pkin(4) * qJD(4);
	t257 = qJD(4) + qJD(5);
	t300 = pkin(5) * t257;
	t238 = -t248 * t300 - t251 * t294;
	t301 = pkin(5) * t248;
	t242 = -pkin(4) * t251 - t301;
	t258 = sin(qJ(2));
	t260 = cos(qJ(2));
	t283 = t258 * qJD(3);
	t295 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8) + qJ(3);
	t305 = t295 * t260;
	t309 = (-t237 * t258 + t305) * qJD(2) + (pkin(7) + sin(pkin(11)) * pkin(3) - t242) * qJD(1) + t260 * t238 + t283;
	t250 = qJ(6) + t254;
	t245 = sin(t250);
	t298 = r_i_i_C(2) * t245;
	t246 = cos(t250);
	t299 = r_i_i_C(1) * t246;
	t267 = t237 - t298 + t299;
	t263 = -t267 * t258 + t305;
	t261 = cos(qJ(1));
	t253 = qJD(6) + t257;
	t275 = t253 * t260 - qJD(1);
	t307 = t261 * t275;
	t297 = r_i_i_C(2) * t246;
	t272 = r_i_i_C(1) * t245 + t297;
	t304 = t272 * t253 - t238;
	t288 = qJD(1) * t260;
	t274 = -t253 + t288;
	t259 = sin(qJ(1));
	t286 = qJD(2) * t258;
	t280 = t259 * t286;
	t303 = t274 * t261 - t280;
	t292 = t253 * t258;
	t284 = qJD(2) * t261;
	t279 = t258 * t284;
	t265 = t274 * t259 + t279;
	t233 = t265 * t245 - t246 * t307;
	t234 = t245 * t307 + t265 * t246;
	t291 = t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
	t269 = t275 * t259;
	t235 = t303 * t245 + t246 * t269;
	t236 = t245 * t269 - t303 * t246;
	t290 = -t235 * r_i_i_C(1) + t236 * r_i_i_C(2);
	t289 = qJD(1) * t259;
	t287 = qJD(1) * t261;
	t285 = qJD(2) * t260;
	t282 = t249 * t300;
	t281 = t253 * t299;
	t278 = t295 * t258;
	t273 = -t257 + t288;
	t271 = t242 * t288 - t238;
	t270 = t249 * (-t257 * t260 + qJD(1));
	t239 = t252 * t294 + t282;
	t266 = qJD(1) * t243 - t239 * t260 - t242 * t286;
	t264 = t239 + (-t237 * t260 - pkin(1) - t278) * qJD(1);
	t262 = qJD(3) * t260 + t304 * t258 + (-t267 * t260 - t278) * qJD(2);
	t241 = t292 * t298;
	t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t309 * t259 + t264 * t261, t262 * t261 - t263 * t289, -t258 * t289 + t260 * t284, -t271 * t259 + t266 * t261 + t291, (t261 * t270 + (t273 * t259 + t279) * t248) * pkin(5) + t291, t291; -t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t264 * t259 + t309 * t261, t262 * t259 + t263 * t287, t258 * t287 + t259 * t285, t266 * t259 + t271 * t261 + t290, (t259 * t270 + (-t273 * t261 + t280) * t248) * pkin(5) + t290, t290; 0, t263 * qJD(2) - t304 * t260 + t283, t286, t241 + (-t239 - t281) * t258 + (t242 - t272) * t285, t241 + (-t281 - t282) * t258 + (-t272 - t301) * t285, -t285 * t297 + t241 + (-t245 * t285 - t246 * t292) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end