% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:12
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
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(7) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0, 0; 0, 0, -t29, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:14
	% EndTime: 2019-10-10 08:54:14
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (87->40), mult. (278->71), div. (0->0), fcn. (217->6), ass. (0->32)
	t201 = sin(qJ(3));
	t204 = cos(qJ(3));
	t223 = pkin(8) + r_i_i_C(3);
	t215 = t223 * t204;
	t229 = -pkin(3) * t201 - qJ(2) + t215;
	t200 = sin(qJ(4));
	t203 = cos(qJ(4));
	t209 = r_i_i_C(1) * t203 - r_i_i_C(2) * t200 + pkin(3);
	t216 = t223 * t201;
	t228 = t209 * t204 + t216;
	t205 = cos(qJ(1));
	t211 = qJD(1) * t201 + qJD(4);
	t226 = t211 * t205;
	t202 = sin(qJ(1));
	t218 = qJD(3) * t205;
	t225 = t211 * t202 - t204 * t218;
	t224 = -pkin(1) - pkin(7);
	t222 = qJD(1) * t202;
	t221 = qJD(1) * t205;
	t220 = qJD(3) * t201;
	t219 = qJD(3) * t204;
	t217 = qJD(4) * t204;
	t212 = -qJD(4) * t201 - qJD(1);
	t210 = r_i_i_C(1) * t200 + r_i_i_C(2) * t203;
	t208 = t212 * t205;
	t207 = t210 * qJD(4);
	t206 = qJD(2) + (pkin(3) * t204 + t216) * qJD(3);
	t199 = t203 * t226 + (t212 * t200 + t203 * t219) * t202;
	t198 = t212 * t203 * t202 + (-t202 * t219 - t226) * t200;
	t197 = t200 * t208 - t225 * t203;
	t196 = t225 * t200 + t203 * t208;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) + t206 * t205 + (t229 * t202 + t224 * t205) * qJD(1), t221, t228 * t221 + (-t210 * t217 + (-t209 * t201 + t215) * qJD(3)) * t202, t198 * r_i_i_C(1) - t199 * r_i_i_C(2), 0, 0; t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t206 * t202 + (t224 * t202 - t229 * t205) * qJD(1), t222, (t209 * t218 + t223 * t222) * t201 + (t209 * t222 + (-t223 * qJD(3) + t207) * t205) * t204, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0; 0, 0, -t228 * qJD(3) + t201 * t207, (t200 * t217 + t203 * t220) * r_i_i_C(2) + (t200 * t220 - t203 * t217) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:14
	% EndTime: 2019-10-10 08:54:14
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (265->53), mult. (428->81), div. (0->0), fcn. (338->8), ass. (0->47)
	t235 = cos(qJ(4));
	t227 = t235 * pkin(4) + pkin(3);
	t233 = sin(qJ(3));
	t236 = cos(qJ(3));
	t266 = r_i_i_C(3) + pkin(9) + pkin(8);
	t251 = t266 * t236;
	t257 = qJD(4) * t235;
	t280 = -pkin(4) * t257 + (-t227 * t233 - qJ(2) + t251) * qJD(1);
	t231 = qJ(4) + qJ(5);
	t228 = sin(t231);
	t229 = cos(t231);
	t276 = -r_i_i_C(1) * t229 + r_i_i_C(2) * t228;
	t242 = t227 - t276;
	t252 = t266 * t233;
	t278 = t242 * t236 + t252;
	t277 = r_i_i_C(1) * t228 + r_i_i_C(2) * t229;
	t275 = t235 * (qJD(4) * t233 + qJD(1));
	t230 = qJD(4) + qJD(5);
	t232 = sin(qJ(4));
	t271 = pkin(4) * t232;
	t256 = qJD(4) * t271;
	t240 = t277 * t230 + t256;
	t234 = sin(qJ(1));
	t262 = qJD(1) * t233;
	t248 = t230 + t262;
	t237 = cos(qJ(1));
	t258 = qJD(3) * t237;
	t253 = t236 * t258;
	t273 = t248 * t234 - t253;
	t254 = qJD(3) * t234 * t236;
	t272 = t248 * t237 + t254;
	t249 = -t230 * t233 - qJD(1);
	t243 = t249 * t237;
	t220 = t273 * t228 + t229 * t243;
	t221 = t228 * t243 - t273 * t229;
	t264 = -t220 * r_i_i_C(1) + t221 * r_i_i_C(2);
	t244 = t249 * t234;
	t222 = -t272 * t228 + t229 * t244;
	t223 = t228 * t244 + t272 * t229;
	t263 = t222 * r_i_i_C(1) - t223 * r_i_i_C(2);
	t261 = qJD(1) * t234;
	t260 = qJD(1) * t237;
	t259 = qJD(3) * t233;
	t246 = -qJD(4) - t262;
	t241 = t276 * t230 * t236 + t277 * t259;
	t239 = -t233 * t256 + qJD(2) + (t227 * t236 + t252) * qJD(3) + (-pkin(1) - pkin(7) - t271) * qJD(1);
	t1 = [t221 * r_i_i_C(1) + t220 * r_i_i_C(2) + t280 * t234 + t239 * t237, t260, t278 * t260 + (-t240 * t236 + (-t242 * t233 + t251) * qJD(3)) * t234, (-t234 * t275 + (t246 * t237 - t254) * t232) * pkin(4) + t263, t263, 0; t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t239 * t234 - t280 * t237, t261, (t242 * t258 + t266 * t261) * t233 + (t242 * t261 + (-t266 * qJD(3) + t240) * t237) * t236, (t237 * t275 + (t246 * t234 + t253) * t232) * pkin(4) + t264, t264, 0; 0, 0, -t278 * qJD(3) + t240 * t233, (t232 * t259 - t236 * t257) * pkin(4) + t241, t241, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:14
	% EndTime: 2019-10-10 08:54:15
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (527->69), mult. (746->98), div. (0->0), fcn. (626->8), ass. (0->55)
	t285 = qJ(4) + qJ(5);
	t282 = sin(t285);
	t284 = qJD(4) + qJD(5);
	t286 = sin(qJ(4));
	t329 = pkin(4) * qJD(4);
	t331 = r_i_i_C(2) + pkin(9) + pkin(8);
	t295 = t331 * qJD(3) + qJD(6) * t282 - t286 * t329;
	t283 = cos(t285);
	t330 = r_i_i_C(3) + qJ(6);
	t309 = t330 * t283;
	t333 = r_i_i_C(1) + pkin(5);
	t341 = (-t333 * t282 + t309) * t284 + t295;
	t289 = cos(qJ(4));
	t281 = pkin(4) * t289 + pkin(3);
	t310 = t330 * t282;
	t301 = -t333 * t283 - t310;
	t340 = -t281 + t301;
	t287 = sin(qJ(3));
	t290 = cos(qJ(3));
	t316 = t289 * t329;
	t317 = t283 * qJD(6);
	t338 = (-t281 * t287 + t290 * t331 - qJ(2)) * qJD(1) - t316 + t317;
	t337 = t289 * (qJD(4) * t287 + qJD(1));
	t332 = pkin(4) * t286;
	t328 = t284 * t287;
	t327 = t284 * t290;
	t288 = sin(qJ(1));
	t326 = t288 * t282;
	t325 = t288 * t283;
	t291 = cos(qJ(1));
	t324 = t291 * t282;
	t323 = qJD(1) * t287;
	t322 = qJD(1) * t288;
	t321 = qJD(1) * t291;
	t320 = qJD(3) * t287;
	t319 = qJD(3) * t290;
	t318 = qJD(3) * t291;
	t315 = t291 * t287 * t283;
	t313 = t282 * t320;
	t314 = t290 * t317 + t313 * t333;
	t312 = t288 * t319;
	t311 = t290 * t318;
	t307 = t284 + t323;
	t305 = -qJD(4) - t323;
	t302 = t287 * t325 + t324;
	t266 = -t282 * t311 - t283 * t321 - t284 * t315 + t307 * t326;
	t267 = -t283 * t311 + (t287 * t324 + t325) * t284 + t302 * qJD(1);
	t300 = -(t315 - t326) * qJD(6) + t330 * t267 - t333 * t266;
	t297 = t291 * t307 + t312;
	t268 = (qJD(1) + t328) * t325 + t297 * t282;
	t269 = -t282 * t322 + t283 * t297 - t326 * t328;
	t299 = qJD(6) * t302 - t268 * t333 + t269 * t330;
	t296 = qJD(3) * t340;
	t293 = t281 * t319 + qJD(2) + (-pkin(1) - pkin(7) - t332) * qJD(1) + t295 * t287;
	t1 = [-t330 * t266 - t333 * t267 + t288 * t338 + t293 * t291, t321, (t287 * t331 - t290 * t340) * t321 + (t287 * t296 + t290 * t341) * t288, (-t288 * t337 + (t291 * t305 - t312) * t286) * pkin(4) + t299, t299, t268; t330 * t268 + t333 * t269 + t293 * t288 - t291 * t338, t322, (-t318 * t340 + t322 * t331) * t287 + (-t291 * t341 - t322 * t340) * t290, (t291 * t337 + (t288 * t305 + t311) * t286) * pkin(4) + t300, t300, t266; 0, 0, -t287 * t341 + t290 * t296, (-t309 + t332) * t320 + (t284 * t301 - t316) * t290 + t314, -t310 * t327 + (-t320 * t330 - t327 * t333) * t283 + t314, t283 * t327 - t313;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end