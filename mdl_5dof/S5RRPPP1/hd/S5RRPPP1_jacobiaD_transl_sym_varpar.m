% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPP1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.14s
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
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:43
	% EndTime: 2019-12-29 18:08:44
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (81->44), mult. (292->79), div. (0->0), fcn. (251->8), ass. (0->33)
	t205 = sin(pkin(8));
	t207 = cos(pkin(8));
	t209 = sin(qJ(2));
	t208 = cos(pkin(5));
	t211 = cos(qJ(2));
	t232 = t208 * t211;
	t216 = -(t205 * t232 + t207 * t209) * r_i_i_C(1) - (-t205 * t209 + t207 * t232) * r_i_i_C(2) - t209 * pkin(2);
	t206 = sin(pkin(5));
	t227 = t206 * qJD(3);
	t222 = t209 * t227;
	t234 = t206 * t211;
	t240 = (qJ(3) * t234 + t216) * qJD(2) + t222;
	t219 = r_i_i_C(1) * t205 + r_i_i_C(2) * t207;
	t238 = t208 * qJ(3) + t219 * t206 + pkin(7);
	t236 = r_i_i_C(3) + qJ(3);
	t235 = t206 * t209;
	t233 = t208 * t209;
	t210 = sin(qJ(1));
	t231 = t210 * t208;
	t230 = qJD(1) * t210;
	t229 = qJD(2) * t206;
	t212 = cos(qJ(1));
	t228 = qJD(2) * t212;
	t226 = t212 * qJD(3);
	t224 = t208 * t230;
	t223 = t211 * t229;
	t221 = t236 * t206;
	t215 = (t205 * t233 - t207 * t211) * r_i_i_C(1) + (t205 * t211 + t207 * t233) * r_i_i_C(2) - t211 * pkin(2);
	t214 = -qJ(3) * t235 - pkin(1) + t215;
	t213 = t236 * t234 + t216;
	t204 = t210 * t223 + (t212 * t235 + t231) * qJD(1);
	t203 = t212 * t223 + (t212 * t208 - t210 * t235) * qJD(1);
	t1 = [t208 * t226 - t204 * r_i_i_C(3) - t240 * t210 + (-t238 * t210 + t214 * t212) * qJD(1), ((t207 * r_i_i_C(1) - t205 * r_i_i_C(2) + pkin(2)) * t230 + (t219 * t208 - t221) * t228) * t209 + ((t205 * t224 - t207 * t228) * r_i_i_C(1) + (t205 * t228 + t207 * t224) * r_i_i_C(2) - pkin(2) * t228 + (-t236 * t230 + t226) * t206) * t211, t203, 0, 0; qJD(3) * t231 + t203 * r_i_i_C(3) + t240 * t212 + (t214 * t210 + t238 * t212) * qJD(1), t213 * t212 * qJD(1) + (t211 * t227 + (-t209 * t221 + t215) * qJD(2)) * t210, t204, 0, 0; 0, qJD(2) * t213 + t222, t209 * t229, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:44
	% EndTime: 2019-12-29 18:08:44
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (142->59), mult. (507->104), div. (0->0), fcn. (463->8), ass. (0->48)
	t246 = sin(qJ(2));
	t248 = cos(qJ(2));
	t243 = sin(pkin(5));
	t266 = t243 * (r_i_i_C(1) + qJ(3));
	t257 = -t246 * pkin(2) + t248 * t266;
	t289 = t257 * qJD(1);
	t245 = cos(pkin(5));
	t273 = qJD(2) * t243;
	t268 = t248 * t273;
	t269 = t243 * qJD(3);
	t244 = cos(pkin(8));
	t270 = qJD(4) * t244;
	t242 = sin(pkin(8));
	t276 = t248 * t242;
	t288 = (qJD(1) * t245 + t268) * qJ(3) + (-qJD(2) * pkin(2) + t245 * t270 + t269) * t246 + qJD(1) * pkin(7) + qJD(4) * t276;
	t287 = r_i_i_C(2) - pkin(3);
	t285 = t248 * pkin(2);
	t283 = r_i_i_C(3) + qJ(4);
	t282 = t242 * t246;
	t247 = sin(qJ(1));
	t281 = t243 * t247;
	t249 = cos(qJ(1));
	t280 = t243 * t249;
	t279 = t244 * t248;
	t278 = t245 * t246;
	t277 = t247 * t245;
	t275 = t248 * t249;
	t274 = t249 * t245;
	t272 = qJD(2) * t247;
	t271 = qJD(2) * t249;
	t267 = t248 * t271;
	t264 = t245 * t279 - t282;
	t263 = t244 * t278 + t276;
	t262 = t244 * t246 + t245 * t276;
	t261 = t242 * t278 - t279;
	t260 = t246 * t274 - t281;
	t259 = t246 * t277 + t280;
	t256 = qJD(1) * t264;
	t255 = qJD(1) * t262;
	t254 = t263 * qJD(2);
	t253 = t262 * qJD(2);
	t251 = -t243 * t270 + t245 * qJD(3) + (-t246 * t243 * qJ(3) - pkin(1) - t285) * qJD(1);
	t250 = t264 * qJD(4) + t248 * t269 + (-t246 * t266 - t285) * qJD(2);
	t238 = t247 * t268 + (t246 * t280 + t277) * qJD(1);
	t237 = t243 * t267 + (-t246 * t281 + t274) * qJD(1);
	t231 = t264 * t272 + (t242 * t275 + t260 * t244) * qJD(1);
	t229 = t271 * t282 - t245 * t244 * t267 + (t259 * t244 + t247 * t276) * qJD(1);
	t1 = [-t238 * r_i_i_C(1) + t287 * (-t262 * t272 + (-t260 * t242 + t244 * t275) * qJD(1)) - t283 * t231 + t251 * t249 - t288 * t247, -t287 * (t247 * t255 + t261 * t271) + t283 * (-t247 * t256 - t263 * t271) - t247 * t289 + t250 * t249, t237, -t229, 0; t237 * r_i_i_C(1) + t287 * (t249 * t253 + (-t259 * t242 + t247 * t279) * qJD(1)) - t283 * t229 + t251 * t247 + t288 * t249, t287 * (t249 * t255 - t261 * t272) - t283 * (t247 * t254 - t249 * t256) + t249 * t289 + t250 * t247, t238, t231, 0; 0, t263 * qJD(4) + t246 * t269 + t287 * t253 + (t264 * t283 + t257) * qJD(2), t246 * t273, t254, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:44
	% EndTime: 2019-12-29 18:08:45
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (194->70), mult. (687->118), div. (0->0), fcn. (638->8), ass. (0->55)
	t271 = sin(pkin(8));
	t272 = sin(pkin(5));
	t273 = cos(pkin(8));
	t274 = cos(pkin(5));
	t275 = sin(qJ(2));
	t277 = cos(qJ(2));
	t292 = qJD(4) * t273 - qJD(5) * t271;
	t302 = t272 * qJD(3);
	t304 = qJD(2) * t277;
	t326 = (-qJD(2) * pkin(2) + t292 * t274 + t302) * t275 + (qJD(1) * t274 + t272 * t304) * qJ(3) + (t271 * qJD(4) + t273 * qJD(5)) * t277 + qJD(1) * pkin(7);
	t323 = -r_i_i_C(1) - pkin(4);
	t296 = t272 * (qJ(3) - t323);
	t281 = -t275 * pkin(2) + t277 * t296;
	t321 = t277 * pkin(2);
	t320 = r_i_i_C(2) + qJ(4);
	t318 = t271 * t275;
	t317 = t272 * t275;
	t276 = sin(qJ(1));
	t316 = t272 * t276;
	t278 = cos(qJ(1));
	t315 = t272 * t278;
	t314 = t273 * t275;
	t313 = t274 * t275;
	t312 = t276 * t274;
	t311 = t277 * t271;
	t310 = t277 * t273;
	t309 = t277 * t278;
	t308 = t278 * t274;
	t307 = qJD(1) * t276;
	t306 = qJD(1) * t278;
	t305 = qJD(2) * t276;
	t303 = qJD(2) * t278;
	t301 = pkin(3) + r_i_i_C(3) + qJ(5);
	t299 = t271 * t313;
	t298 = t276 * t304;
	t297 = t277 * t303;
	t295 = qJD(1) * t299;
	t294 = t271 * t274 * t305;
	t290 = t274 * t310 - t318;
	t289 = t273 * t313 + t311;
	t288 = t274 * t311 + t314;
	t287 = t299 - t310;
	t285 = qJD(1) * t290;
	t284 = qJD(1) * t288;
	t283 = t289 * qJD(2);
	t282 = t288 * qJD(2);
	t280 = t274 * qJD(3) - t292 * t272 + (-qJ(3) * t317 - pkin(1) - t321) * qJD(1);
	t279 = -t288 * qJD(5) + t290 * qJD(4) + t277 * t302 + (-t275 * t296 - t321) * qJD(2);
	t261 = t272 * t298 + (t275 * t315 + t312) * qJD(1);
	t260 = t272 * t297 + (-t275 * t316 + t308) * qJD(1);
	t255 = -t277 * t294 - t278 * t295 - t305 * t314 + (t271 * t316 + t273 * t309) * qJD(1);
	t254 = t290 * t305 + (t271 * t309 + (t275 * t308 - t316) * t273) * qJD(1);
	t253 = -t272 * t271 * t306 - t276 * t295 + t278 * t282 + t307 * t310;
	t252 = t303 * t318 - t274 * t273 * t297 + (t276 * t311 + (t275 * t312 + t315) * t273) * qJD(1);
	t1 = [-t320 * t254 - t301 * t255 + t323 * t261 - t326 * t276 + t280 * t278, t320 * (-t276 * t285 - t289 * t303) + t301 * (t276 * t284 + t287 * t303) - t281 * t307 + t279 * t278, t260, -t252, -t253; -t320 * t252 - t301 * t253 - t323 * t260 + t280 * t276 + t326 * t278, -t320 * (t276 * t283 - t278 * t285) - t301 * (t273 * t298 - t275 * t294 + t278 * t284) + t281 * t306 + t279 * t276, t261, t254, t255; 0, t289 * qJD(4) - t287 * qJD(5) + t275 * t302 - t301 * t282 + (t290 * t320 + t281) * qJD(2), qJD(2) * t317, t283, -t287 * qJD(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end