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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:38
	% EndTime: 2019-12-31 19:25:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.07s
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:40
	% EndTime: 2019-12-31 19:25:40
	% DurationCPUTime: 0.28s
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
	t1 = [t208 * t226 - t204 * r_i_i_C(3) - t240 * t210 + (-t238 * t210 + t214 * t212) * qJD(1), ((t207 * r_i_i_C(1) - t205 * r_i_i_C(2) + pkin(2)) * t230 + (t219 * t208 - t221) * t228) * t209 + ((t205 * t224 - t207 * t228) * r_i_i_C(1) + (t205 * t228 + t207 * t224) * r_i_i_C(2) - pkin(2) * t228 + (-t236 * t230 + t226) * t206) * t211, t203, 0, 0; qJD(3) * t231 + t203 * r_i_i_C(3) + t240 * t212 + (t214 * t210 + t238 * t212) * qJD(1), t213 * t212 * qJD(1) + (t211 * t227 + (-t209 * t221 + t215) * qJD(2)) * t210, t204, 0, 0; 0, t213 * qJD(2) + t222, t209 * t229, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:40
	% EndTime: 2019-12-31 19:25:40
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (142->59), mult. (507->105), div. (0->0), fcn. (463->8), ass. (0->47)
	t243 = sin(qJ(2));
	t245 = cos(qJ(2));
	t240 = sin(pkin(5));
	t263 = t240 * (r_i_i_C(1) + qJ(3));
	t254 = -t243 * pkin(2) + t245 * t263;
	t285 = t254 * qJD(1);
	t242 = cos(pkin(5));
	t270 = qJD(2) * t240;
	t265 = t245 * t270;
	t266 = t240 * qJD(3);
	t241 = cos(pkin(8));
	t267 = qJD(4) * t241;
	t239 = sin(pkin(8));
	t273 = t245 * t239;
	t284 = (qJD(1) * t242 + t265) * qJ(3) + (-qJD(2) * pkin(2) + t242 * t267 + t266) * t243 + qJD(1) * pkin(7) + qJD(4) * t273;
	t283 = -r_i_i_C(2) + pkin(3);
	t281 = t245 * pkin(2);
	t279 = r_i_i_C(3) + qJ(4);
	t278 = t239 * t243;
	t246 = cos(qJ(1));
	t277 = t240 * t246;
	t276 = t241 * t245;
	t275 = t242 * t243;
	t244 = sin(qJ(1));
	t274 = t243 * t244;
	t272 = t245 * t246;
	t271 = t246 * t242;
	t269 = qJD(2) * t244;
	t268 = qJD(2) * t246;
	t264 = t245 * t268;
	t261 = t242 * t276 - t278;
	t260 = t241 * t275 + t273;
	t259 = t241 * t243 + t242 * t273;
	t258 = t239 * t275 - t276;
	t257 = -t240 * t244 + t243 * t271;
	t256 = t242 * t274 + t277;
	t253 = qJD(1) * t261;
	t252 = qJD(1) * t259;
	t251 = t260 * qJD(2);
	t250 = t259 * qJD(2);
	t248 = -t240 * t267 + t242 * qJD(3) + (-qJ(3) * t240 * t243 - pkin(1) - t281) * qJD(1);
	t247 = t261 * qJD(4) + t245 * t266 + (-t243 * t263 - t281) * qJD(2);
	t235 = t244 * t265 + (t242 * t244 + t243 * t277) * qJD(1);
	t234 = t240 * t264 + (-t240 * t274 + t271) * qJD(1);
	t228 = t261 * t269 + (t239 * t272 + t257 * t241) * qJD(1);
	t226 = t268 * t278 - t242 * t241 * t264 + (t256 * t241 + t244 * t273) * qJD(1);
	t1 = [-t235 * r_i_i_C(1) - t283 * (-t259 * t269 + (-t257 * t239 + t241 * t272) * qJD(1)) - t279 * t228 + t248 * t246 - t284 * t244, t283 * (t244 * t252 + t258 * t268) + t279 * (-t244 * t253 - t260 * t268) - t244 * t285 + t247 * t246, t234, -t226, 0; t234 * r_i_i_C(1) - t283 * (t246 * t250 + (-t256 * t239 + t244 * t276) * qJD(1)) - t279 * t226 + t248 * t244 + t284 * t246, -t283 * (t246 * t252 - t258 * t269) - t279 * (t244 * t251 - t246 * t253) + t246 * t285 + t247 * t244, t235, t228, 0; 0, t260 * qJD(4) + t243 * t266 - t283 * t250 + (t261 * t279 + t254) * qJD(2), t243 * t270, t251, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:40
	% EndTime: 2019-12-31 19:25:40
	% DurationCPUTime: 0.33s
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
	t323 = pkin(4) + r_i_i_C(1);
	t296 = t272 * (qJ(3) + t323);
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
	t300 = -r_i_i_C(3) - qJ(5) - pkin(3);
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
	t1 = [-t320 * t254 + t300 * t255 - t323 * t261 - t326 * t276 + t280 * t278, t320 * (-t276 * t285 - t289 * t303) - t300 * (t276 * t284 + t287 * t303) - t281 * t307 + t279 * t278, t260, -t252, -t253; -t320 * t252 + t300 * t253 + t323 * t260 + t280 * t276 + t326 * t278, -t320 * (t276 * t283 - t278 * t285) + t300 * (t273 * t298 - t275 * t294 + t278 * t284) + t281 * t306 + t279 * t276, t261, t254, t255; 0, t289 * qJD(4) - t287 * qJD(5) + t275 * t302 + t300 * t282 + (t290 * t320 + t281) * qJD(2), qJD(2) * t317, t283, -t287 * qJD(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end