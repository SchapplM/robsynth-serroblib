% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (153->22), mult. (174->40), div. (0->0), fcn. (146->8), ass. (0->23)
	t170 = sin(pkin(5));
	t185 = t170 * (pkin(8) + r_i_i_C(3));
	t183 = pkin(1) * qJD(1);
	t171 = cos(pkin(5));
	t172 = sin(qJ(3));
	t182 = t171 * t172;
	t173 = cos(qJ(3));
	t181 = t171 * t173;
	t169 = qJ(1) + qJ(2);
	t166 = sin(t169);
	t167 = cos(t169);
	t179 = t166 * t172 - t167 * t181;
	t178 = t166 * t173 + t167 * t182;
	t177 = t166 * t181 + t167 * t172;
	t176 = t166 * t182 - t167 * t173;
	t168 = qJD(1) + qJD(2);
	t160 = t176 * qJD(3) + t179 * t168;
	t161 = t177 * qJD(3) + t178 * t168;
	t175 = -t161 * r_i_i_C(1) + t160 * r_i_i_C(2) + (-pkin(2) * t166 + t167 * t185) * t168;
	t162 = t178 * qJD(3) + t177 * t168;
	t163 = t179 * qJD(3) + t176 * t168;
	t174 = t162 * r_i_i_C(2) + t163 * r_i_i_C(1) + (-pkin(2) * t167 - t166 * t185) * t168;
	t1 = [-cos(qJ(1)) * t183 + t174, t174, t160 * r_i_i_C(1) + t161 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t183 + t175, t175, -t162 * r_i_i_C(1) + t163 * r_i_i_C(2), 0, 0; 0, 0, (-r_i_i_C(1) * t172 - r_i_i_C(2) * t173) * t170 * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (509->48), mult. (351->72), div. (0->0), fcn. (244->14), ass. (0->47)
	t255 = qJD(3) + qJD(4);
	t285 = -t255 / 0.2e1;
	t257 = qJ(3) + qJ(4);
	t250 = pkin(5) - t257;
	t246 = sin(t250);
	t288 = t246 * t285;
	t249 = pkin(5) + t257;
	t283 = sin(t249);
	t235 = t283 / 0.2e1 - t246 / 0.2e1;
	t258 = qJ(1) + qJ(2);
	t252 = sin(t258);
	t254 = cos(t258);
	t256 = qJD(1) + qJD(2);
	t284 = cos(t250);
	t286 = cos(t249) / 0.2e1;
	t236 = t284 / 0.2e1 + t286;
	t253 = cos(t257);
	t270 = t236 * t255 + t253 * t256;
	t251 = sin(t257);
	t280 = t251 * t255;
	t287 = t270 * t252 + (t235 * t256 + t280) * t254;
	t282 = pkin(1) * qJD(1);
	t281 = pkin(3) * qJD(3);
	t279 = t252 * t256;
	t260 = cos(pkin(5));
	t261 = sin(qJ(3));
	t277 = t260 * t261;
	t262 = cos(qJ(3));
	t276 = t260 * t262;
	t267 = t236 * t256 + t253 * t255;
	t240 = t283 * t285;
	t269 = t251 * t256 - t240 + t288;
	t228 = t269 * t252 - t267 * t254;
	t275 = t228 * r_i_i_C(1) + t287 * r_i_i_C(2);
	t229 = t267 * t252 + t269 * t254;
	t265 = t235 * t279 + t252 * t280 - t270 * t254;
	t274 = -t229 * r_i_i_C(1) + t265 * r_i_i_C(2);
	t273 = (t255 * t286 + t284 * t285) * r_i_i_C(1) + (t240 + t288) * r_i_i_C(2);
	t259 = sin(pkin(5));
	t272 = r_i_i_C(3) * t256 * t259;
	t271 = t261 * t281;
	t266 = -t252 * t276 - t254 * t261;
	t237 = pkin(3) * t277 - t259 * (pkin(8) + pkin(9));
	t248 = t262 * pkin(3) + pkin(2);
	t264 = t229 * r_i_i_C(2) + t265 * r_i_i_C(1) + t237 * t279 + (-t248 * t256 - t276 * t281) * t254 + (-t272 + t271) * t252;
	t263 = t228 * r_i_i_C(2) - t287 * r_i_i_C(1) + t254 * t272 + (-t237 * t254 - t248 * t252) * t256 + t266 * t281;
	t1 = [-cos(qJ(1)) * t282 + t264, t264, ((t252 * t261 - t254 * t276) * t256 + (t252 * t277 - t254 * t262) * qJD(3)) * pkin(3) + t275, t275, 0; -sin(qJ(1)) * t282 + t263, t263, (t266 * t256 + (-t252 * t262 - t254 * t277) * qJD(3)) * pkin(3) + t274, t274, 0; 0, 0, -t259 * t271 + t273, t273, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (938->71), mult. (628->104), div. (0->0), fcn. (454->18), ass. (0->65)
	t278 = qJD(3) + qJD(4);
	t272 = qJD(5) + t278;
	t320 = -t272 / 0.2e1;
	t279 = qJD(1) + qJD(2);
	t282 = sin(pkin(5));
	t285 = sin(qJ(3));
	t305 = qJD(3) * t285;
	t280 = qJ(3) + qJ(4);
	t273 = sin(t280);
	t317 = pkin(4) * t273;
	t324 = r_i_i_C(3) * t279 * t282 - pkin(3) * t305 - t278 * t317;
	t277 = qJ(5) + t280;
	t268 = pkin(5) - t277;
	t265 = sin(t268);
	t323 = t265 * t320;
	t267 = pkin(5) + t277;
	t318 = sin(t267);
	t252 = t318 / 0.2e1 - t265 / 0.2e1;
	t281 = qJ(1) + qJ(2);
	t274 = sin(t281);
	t276 = cos(t281);
	t319 = cos(t268);
	t321 = cos(t267) / 0.2e1;
	t253 = t319 / 0.2e1 + t321;
	t270 = cos(t277);
	t300 = t253 * t272 + t270 * t279;
	t269 = sin(t277);
	t315 = t269 * t272;
	t322 = t300 * t274 + (t252 * t279 + t315) * t276;
	t316 = pkin(1) * qJD(1);
	t313 = t274 * t279;
	t275 = cos(t280);
	t312 = t275 * t278;
	t311 = t276 * t279;
	t284 = sin(qJ(4));
	t310 = t284 * t285;
	t287 = cos(qJ(3));
	t309 = t284 * t287;
	t296 = t253 * t279 + t270 * t272;
	t257 = t318 * t320;
	t299 = t269 * t279 - t257 + t323;
	t238 = t299 * t274 - t296 * t276;
	t308 = t238 * r_i_i_C(1) + t322 * r_i_i_C(2);
	t239 = t296 * t274 + t299 * t276;
	t292 = t252 * t313 + t274 * t315 - t300 * t276;
	t307 = -t239 * r_i_i_C(1) + t292 * r_i_i_C(2);
	t306 = (t272 * t321 + t319 * t320) * r_i_i_C(1) + (t257 + t323) * r_i_i_C(2);
	t304 = qJD(3) * t287;
	t283 = cos(pkin(5));
	t286 = cos(qJ(4));
	t271 = t286 * pkin(4) + pkin(3);
	t295 = -t285 * t286 - t309;
	t293 = t295 * qJD(4);
	t288 = -t271 * t305 + (-t284 * t304 + t293) * pkin(4);
	t302 = (-t285 * pkin(3) - t317) * t279 + t288 * t283;
	t294 = t286 * t287 - t310;
	t301 = -(t287 * pkin(3) + pkin(4) * t275 + pkin(2)) * t279 - (t271 * t304 + (t294 * qJD(4) - t284 * t305) * pkin(4)) * t283;
	t298 = -(-pkin(4) * t310 + t271 * t287) * t283 * t279 - pkin(3) * t304 - pkin(4) * t312;
	t244 = -t282 * (pkin(8) + pkin(9) + pkin(10)) + (pkin(4) * t309 + t285 * t271) * t283;
	t291 = -t322 * r_i_i_C(1) + t238 * r_i_i_C(2) - t244 * t311 + t301 * t274 + t324 * t276;
	t290 = t292 * r_i_i_C(1) + t239 * r_i_i_C(2) + t244 * t313 - t324 * t274 + t301 * t276;
	t289 = pkin(4) * (t295 * qJD(3) + t293);
	t251 = t294 * t283 * pkin(4);
	t242 = t283 * t289;
	t1 = [-cos(qJ(1)) * t316 + t290, t290, -t302 * t274 + t298 * t276 + t308, -t251 * t311 - t274 * t242 + (t273 * t313 - t276 * t312) * pkin(4) + t308, t308; -sin(qJ(1)) * t316 + t291, t291, t298 * t274 + t302 * t276 + t307, -t251 * t313 + t276 * t242 + (-t273 * t311 - t274 * t312) * pkin(4) + t307, t307; 0, 0, t288 * t282 + t306, t282 * t289 + t306, t306;];
	JaD_transl = t1;
end