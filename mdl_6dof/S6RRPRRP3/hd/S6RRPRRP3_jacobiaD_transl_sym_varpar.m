% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP3
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
% Datum: 2019-10-10 10:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 10:33:37
	% EndTime: 2019-10-10 10:33:38
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
	t212 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:33:38
	% EndTime: 2019-10-10 10:33:38
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (374->56), mult. (446->81), div. (0->0), fcn. (351->10), ass. (0->52)
	t255 = cos(qJ(4));
	t242 = t255 * pkin(4) + pkin(3);
	t249 = qJ(2) + pkin(10);
	t244 = sin(t249);
	t245 = cos(t249);
	t289 = r_i_i_C(3) + pkin(9) + pkin(8);
	t264 = t289 * t245 - sin(qJ(2)) * pkin(2);
	t270 = qJD(4) * t245 - qJD(1);
	t252 = sin(qJ(4));
	t293 = pkin(4) * t252;
	t302 = (-t242 * t244 + t264) * qJD(2) - t270 * t293 - qJD(1) * (-qJ(3) - pkin(7));
	t250 = qJ(4) + qJ(5);
	t247 = cos(t250);
	t246 = sin(t250);
	t292 = r_i_i_C(2) * t246;
	t265 = r_i_i_C(1) * t247 + t242 - t292;
	t261 = -t265 * t244 + t264;
	t257 = cos(qJ(1));
	t248 = qJD(4) + qJD(5);
	t272 = t245 * t248 - qJD(1);
	t300 = t257 * t272;
	t298 = -t289 * t244 - cos(qJ(2)) * pkin(2);
	t291 = r_i_i_C(2) * t247;
	t268 = r_i_i_C(1) * t246 + t291;
	t288 = pkin(4) * qJD(4);
	t297 = t268 * t248 + t252 * t288;
	t283 = qJD(1) * t245;
	t271 = -t248 + t283;
	t254 = sin(qJ(1));
	t279 = qJD(2) * t244;
	t275 = t254 * t279;
	t296 = t271 * t257 - t275;
	t286 = t247 * t248;
	t274 = t257 * t279;
	t262 = t271 * t254 + t274;
	t237 = t262 * t246 - t247 * t300;
	t238 = t246 * t300 + t262 * t247;
	t285 = t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
	t267 = t272 * t254;
	t239 = t296 * t246 + t247 * t267;
	t240 = t246 * t267 - t296 * t247;
	t284 = -t239 * r_i_i_C(1) + t240 * r_i_i_C(2);
	t281 = qJD(1) * t254;
	t280 = qJD(1) * t257;
	t278 = qJD(2) * t245;
	t276 = t255 * t288;
	t269 = -qJD(4) + t283;
	t266 = t270 * t255;
	t260 = t276 + qJD(3) + (-t242 * t245 - pkin(1) + t298) * qJD(1);
	t259 = t297 * t244 + (-t265 * t245 + t298) * qJD(2);
	t241 = t244 * t248 * t292;
	t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t302 * t254 + t260 * t257, t259 * t257 - t261 * t281, t280, (-t257 * t266 + (t269 * t254 + t274) * t252) * pkin(4) + t285, t285, 0; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t260 * t254 + t302 * t257, t259 * t254 + t261 * t280, t281, (-t254 * t266 + (-t269 * t257 + t275) * t252) * pkin(4) + t284, t284, 0; 0, t261 * qJD(2) - t297 * t245, 0, t241 + (-r_i_i_C(1) * t286 - t276) * t244 + (-t268 - t293) * t278, -t278 * t291 + t241 + (-t244 * t286 - t246 * t278) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:38
	% EndTime: 2019-10-10 10:33:38
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (509->69), mult. (548->96), div. (0->0), fcn. (429->10), ass. (0->56)
	t260 = sin(qJ(4));
	t258 = qJ(4) + qJ(5);
	t252 = sin(t258);
	t256 = qJD(4) + qJD(5);
	t296 = t252 * t256;
	t298 = pkin(4) * qJD(4);
	t242 = -pkin(5) * t296 - t260 * t298;
	t253 = cos(t258);
	t263 = cos(qJ(4));
	t247 = pkin(4) * t263 + pkin(5) * t253;
	t245 = pkin(3) + t247;
	t246 = pkin(4) * t260 + pkin(5) * t252;
	t257 = qJ(2) + pkin(10);
	t250 = sin(t257);
	t251 = cos(t257);
	t299 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
	t273 = t299 * t251 - sin(qJ(2)) * pkin(2);
	t284 = t250 * qJD(6);
	t311 = (-t245 * t250 + t273) * qJD(2) + (t246 + qJ(3) + pkin(7)) * qJD(1) + t251 * t242 + t284;
	t265 = cos(qJ(1));
	t291 = qJD(1) * t251;
	t280 = -t256 + t291;
	t262 = sin(qJ(1));
	t286 = qJD(2) * t262;
	t305 = -t250 * t286 + t265 * t280;
	t310 = t305 * t252;
	t274 = r_i_i_C(1) * t253 - r_i_i_C(2) * t252 + t245;
	t267 = -t274 * t250 + t273;
	t307 = -t299 * t250 - cos(qJ(2)) * pkin(2);
	t301 = r_i_i_C(2) * t253;
	t279 = r_i_i_C(1) * t252 + t301;
	t306 = t256 * t279 - t242;
	t303 = -pkin(5) - r_i_i_C(1);
	t295 = t253 * t256;
	t285 = qJD(2) * t265;
	t269 = t250 * t285 + t262 * t280;
	t281 = t251 * t256 - qJD(1);
	t277 = t253 * t281;
	t238 = t252 * t269 - t265 * t277;
	t239 = t252 * t265 * t281 + t253 * t269;
	t294 = r_i_i_C(1) * t238 + r_i_i_C(2) * t239;
	t276 = t281 * t262;
	t240 = t253 * t276 + t310;
	t241 = t252 * t276 - t253 * t305;
	t293 = -r_i_i_C(1) * t240 + r_i_i_C(2) * t241;
	t290 = qJD(1) * t262;
	t289 = qJD(1) * t265;
	t288 = qJD(2) * t250;
	t287 = qJD(2) * t251;
	t278 = t246 * t291 + t242;
	t243 = pkin(5) * t295 + t263 * t298;
	t271 = qJD(1) * t247 - t243 * t251 + t246 * t288;
	t268 = qJD(3) + t243 + (-t245 * t251 - pkin(1) + t307) * qJD(1);
	t266 = qJD(6) * t251 + t306 * t250 + (-t251 * t274 + t307) * qJD(2);
	t244 = t250 * r_i_i_C(2) * t296;
	t1 = [t241 * r_i_i_C(1) + t240 * r_i_i_C(2) - t262 * t311 + t268 * t265, t266 * t265 - t267 * t290, t289, t262 * t278 + t265 * t271 + t294, pkin(5) * t238 + t294, -t250 * t290 + t251 * t285; -t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t268 * t262 + t265 * t311, t262 * t266 + t267 * t289, t290, t262 * t271 - t265 * t278 + t293, (-t262 * t277 - t310) * pkin(5) + t293, t250 * t289 + t251 * t286; 0, qJD(2) * t267 - t251 * t306 + t284, 0, t244 + (-r_i_i_C(1) * t295 - t243) * t250 + (-t246 - t279) * t287, t244 + t303 * t250 * t295 + (t252 * t303 - t301) * t287, t288;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end