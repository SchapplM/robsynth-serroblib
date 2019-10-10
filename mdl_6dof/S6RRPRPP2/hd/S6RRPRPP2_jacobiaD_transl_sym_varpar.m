% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
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
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
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
	% StartTime: 2019-10-10 09:57:25
	% EndTime: 2019-10-10 09:57:25
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(9);
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
	% StartTime: 2019-10-10 09:57:26
	% EndTime: 2019-10-10 09:57:26
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
	t212 = qJ(2) + pkin(9);
	t210 = sin(t212);
	t211 = cos(t212);
	t243 = pkin(8) + r_i_i_C(3);
	t225 = t243 * t211 - sin(qJ(2)) * pkin(2);
	t246 = -pkin(3) * t210 + t225;
	t214 = sin(qJ(4));
	t217 = cos(qJ(4));
	t227 = r_i_i_C(1) * t217 - r_i_i_C(2) * t214 + pkin(3);
	t221 = -t227 * t210 + t225;
	t244 = -t243 * t210 - cos(qJ(2)) * pkin(2);
	t219 = cos(qJ(1));
	t239 = t217 * t219;
	t216 = sin(qJ(1));
	t238 = qJD(1) * t216;
	t237 = qJD(1) * t219;
	t236 = qJD(2) * t216;
	t235 = qJD(2) * t217;
	t234 = qJD(2) * t219;
	t233 = qJD(4) * t210;
	t232 = qJD(4) * t211;
	t230 = -qJD(1) + t232;
	t229 = qJD(1) * t211 - qJD(4);
	t228 = r_i_i_C(1) * t214 + r_i_i_C(2) * t217;
	t226 = t230 * t214;
	t223 = -pkin(3) * t211 - pkin(1) + t244;
	t222 = t210 * t234 + t229 * t216;
	t220 = t228 * t233 + (-t227 * t211 + t244) * qJD(2);
	t213 = -qJ(3) - pkin(7);
	t208 = -t229 * t239 + (t210 * t235 + t226) * t216;
	t207 = t230 * t217 * t216 + (-t210 * t236 + t229 * t219) * t214;
	t206 = t222 * t217 + t219 * t226;
	t205 = t222 * t214 - t230 * t239;
	t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t219 * qJD(3) - t246 * t236 + (t213 * t216 + t223 * t219) * qJD(1), t220 * t219 - t221 * t238, t237, t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0, 0; -t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t216 * qJD(3) + t246 * t234 + (-t213 * t219 + t223 * t216) * qJD(1), t220 * t216 + t221 * t237, t238, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2), 0, 0; 0, t221 * qJD(2) - t228 * t232, 0, (-t211 * t235 + t214 * t233) * r_i_i_C(2) + (-qJD(2) * t211 * t214 - t217 * t233) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:27
	% EndTime: 2019-10-10 09:57:27
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (308->56), mult. (536->83), div. (0->0), fcn. (446->8), ass. (0->45)
	t262 = qJ(2) + pkin(9);
	t260 = sin(t262);
	t261 = cos(t262);
	t304 = pkin(8) + r_i_i_C(2);
	t281 = t304 * t261 - sin(qJ(2)) * pkin(2);
	t308 = -pkin(3) * t260 + t281;
	t264 = sin(qJ(4));
	t267 = cos(qJ(4));
	t299 = r_i_i_C(3) + qJ(5);
	t303 = -r_i_i_C(1) - pkin(4);
	t276 = -t299 * t264 + t303 * t267;
	t272 = -pkin(3) + t276;
	t271 = t272 * t260 + t281;
	t306 = -t304 * t260 - cos(qJ(2)) * pkin(2);
	t275 = t303 * t264 + t299 * t267;
	t305 = t275 * qJD(4) + qJD(5) * t264;
	t266 = sin(qJ(1));
	t298 = t266 * t264;
	t297 = t266 * t267;
	t269 = cos(qJ(1));
	t296 = t269 * t264;
	t295 = t269 * t267;
	t294 = qJD(1) * t266;
	t293 = qJD(1) * t269;
	t292 = qJD(2) * t261;
	t291 = qJD(2) * t266;
	t290 = qJD(2) * t269;
	t289 = qJD(4) * t267;
	t288 = qJD(4) * t269;
	t285 = t260 * t291;
	t284 = qJD(4) * t298;
	t283 = t260 * t290;
	t282 = t267 * t288;
	t279 = t261 * t295 + t298;
	t278 = t261 * t298 + t295;
	t277 = -pkin(3) * t261 - pkin(1) + t306;
	t274 = t264 * t288 + t267 * t294;
	t273 = t264 * t293 + t266 * t289;
	t270 = -t305 * t260 + (t272 * t261 + t306) * qJD(2);
	t263 = -qJ(3) - pkin(7);
	t248 = t279 * qJD(1) - t261 * t284 - t267 * t285 - t282;
	t247 = t273 * t261 - t264 * t285 - t274;
	t246 = t274 * t261 + t267 * t283 - t273;
	t245 = t278 * qJD(1) - t261 * t282 + t264 * t283 - t284;
	t1 = [-t278 * qJD(5) + t269 * qJD(3) + t303 * t248 - t299 * t247 - t308 * t291 + (t266 * t263 + t277 * t269) * qJD(1), t270 * t269 - t271 * t294, t293, t279 * qJD(5) - t303 * t245 - t299 * t246, -t245, 0; -(-t261 * t296 + t297) * qJD(5) + t266 * qJD(3) + t303 * t246 - t299 * t245 + t308 * t290 + (-t269 * t263 + t277 * t266) * qJD(1), t270 * t266 + t271 * t293, t294, -(-t261 * t297 + t296) * qJD(5) + t299 * t248 + t303 * t247, t247, 0; 0, t271 * qJD(2) + t305 * t261, 0, t275 * t292 + (t276 * qJD(4) + t267 * qJD(5)) * t260, t260 * t289 + t264 * t292, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:57:26
	% EndTime: 2019-10-10 09:57:27
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (396->59), mult. (672->85), div. (0->0), fcn. (559->8), ass. (0->44)
	t223 = qJ(2) + pkin(9);
	t221 = sin(t223);
	t222 = cos(t223);
	t248 = pkin(8) - r_i_i_C(3) - qJ(6);
	t241 = t248 * t222 - sin(qJ(2)) * pkin(2);
	t251 = t221 * qJD(6);
	t225 = sin(qJ(4));
	t252 = qJD(5) * t225;
	t271 = (-pkin(3) * t221 + t241) * qJD(2) - qJD(1) * (-qJ(3) - pkin(7)) + t222 * t252 - t251;
	t228 = cos(qJ(4));
	t249 = pkin(4) + pkin(5) + r_i_i_C(1);
	t262 = r_i_i_C(2) + qJ(5);
	t237 = -t262 * t225 - t249 * t228;
	t235 = -pkin(3) + t237;
	t232 = t235 * t221 + t241;
	t266 = t249 * t225 - t262 * t228;
	t269 = t266 * qJD(4) - t252;
	t267 = -t248 * t221 - cos(qJ(2)) * pkin(2);
	t227 = sin(qJ(1));
	t261 = t227 * t225;
	t230 = cos(qJ(1));
	t260 = t230 * t228;
	t259 = qJD(1) * t227;
	t258 = qJD(1) * t230;
	t257 = qJD(2) * t222;
	t256 = qJD(2) * t227;
	t255 = qJD(2) * t230;
	t254 = qJD(4) * t228;
	t253 = qJD(4) * t230;
	t250 = t228 * qJD(5);
	t247 = t221 * t256;
	t246 = qJD(4) * t261;
	t245 = t221 * t255;
	t244 = t228 * t253;
	t242 = t222 * t260 + t261;
	t239 = t225 * t253 + t228 * t259;
	t238 = t225 * t258 + t227 * t254;
	t233 = -t250 + qJD(3) + (-pkin(3) * t222 - pkin(1) + t267) * qJD(1);
	t231 = -qJD(6) * t222 + t269 * t221 + (t235 * t222 + t267) * qJD(2);
	t209 = t242 * qJD(1) - t222 * t246 - t228 * t247 - t244;
	t208 = t238 * t222 - t225 * t247 - t239;
	t207 = t239 * t222 + t228 * t245 - t238;
	t206 = t225 * t245 - t222 * t244 - t246 + (t222 * t261 + t260) * qJD(1);
	t1 = [-t262 * t208 - t249 * t209 - t271 * t227 + t233 * t230, t231 * t230 - t232 * t259, t258, t242 * qJD(5) + t249 * t206 - t262 * t207, -t206, t221 * t259 - t222 * t255; -t262 * t206 - t249 * t207 + t233 * t227 + t271 * t230, t231 * t227 + t232 * t258, t259, -(-t227 * t222 * t228 + t230 * t225) * qJD(5) + t262 * t209 - t249 * t208, t208, -t221 * t258 - t222 * t256; 0, t232 * qJD(2) - t269 * t222 - t251, 0, -t266 * t257 + (t237 * qJD(4) + t250) * t221, t221 * t254 + t225 * t257, -qJD(2) * t221;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end