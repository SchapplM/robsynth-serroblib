% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:16
	% EndTime: 2019-10-09 22:07:16
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(10));
	t182 = cos(pkin(10));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:16
	% EndTime: 2019-10-09 22:07:16
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (116->40), mult. (277->78), div. (0->0), fcn. (252->10), ass. (0->33)
	t218 = r_i_i_C(3) + qJ(4) + pkin(8);
	t195 = sin(pkin(10));
	t196 = sin(pkin(6));
	t217 = t195 * t196;
	t197 = cos(pkin(10));
	t216 = t196 * t197;
	t200 = sin(qJ(3));
	t215 = t196 * t200;
	t201 = sin(qJ(2));
	t214 = t196 * t201;
	t198 = cos(pkin(6));
	t213 = t198 * t201;
	t203 = cos(qJ(2));
	t212 = t198 * t203;
	t211 = qJD(2) * t201;
	t210 = qJD(2) * t203;
	t209 = t195 * t211;
	t208 = t197 * t210;
	t194 = qJ(3) + pkin(11);
	t192 = sin(t194);
	t193 = cos(t194);
	t202 = cos(qJ(3));
	t207 = -pkin(3) * t202 - r_i_i_C(1) * t193 + r_i_i_C(2) * t192 - pkin(2);
	t186 = t195 * t203 + t197 * t213;
	t206 = t195 * t212 + t197 * t201;
	t205 = pkin(3) * t200 + r_i_i_C(1) * t192 + r_i_i_C(2) * t193;
	t204 = qJD(3) * t205;
	t188 = -t195 * t213 + t197 * t203;
	t184 = -t198 * t209 + t208;
	t183 = t206 * qJD(2);
	t182 = t186 * qJD(2);
	t181 = -t198 * t208 + t209;
	t1 = [0, qJD(4) * t188 - t218 * t183 + t207 * t184 + t206 * t204, t205 * t183 + ((-t188 * t193 - t192 * t217) * r_i_i_C(1) + (t188 * t192 - t193 * t217) * r_i_i_C(2) + (-t188 * t202 - t195 * t215) * pkin(3)) * qJD(3), t184, 0, 0; 0, qJD(4) * t186 - t218 * t181 + t207 * t182 - (-t195 * t201 + t197 * t212) * t204, t205 * t181 + ((-t186 * t193 + t192 * t216) * r_i_i_C(1) + (t186 * t192 + t193 * t216) * r_i_i_C(2) + (-t186 * t202 + t197 * t215) * pkin(3)) * qJD(3), t182, 0, 0; 0, (qJD(4) * t201 - t203 * t204 + (t207 * t201 + t218 * t203) * qJD(2)) * t196, -t205 * t196 * t210 + ((-t192 * t198 - t193 * t214) * r_i_i_C(1) + (t192 * t214 - t193 * t198) * r_i_i_C(2) + (-t198 * t200 - t202 * t214) * pkin(3)) * qJD(3), t196 * t211, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:17
	% EndTime: 2019-10-09 22:07:17
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (276->55), mult. (590->97), div. (0->0), fcn. (568->12), ass. (0->43)
	t277 = sin(pkin(10));
	t280 = cos(pkin(10));
	t286 = cos(qJ(2));
	t281 = cos(pkin(6));
	t284 = sin(qJ(2));
	t301 = t281 * t284;
	t266 = t277 * t286 + t280 * t301;
	t275 = qJ(3) + pkin(11);
	t273 = sin(t275);
	t274 = cos(t275);
	t278 = sin(pkin(6));
	t304 = t278 * t280;
	t308 = -t266 * t274 + t273 * t304;
	t307 = r_i_i_C(3) + qJ(5);
	t305 = t277 * t278;
	t283 = sin(qJ(3));
	t303 = t278 * t283;
	t302 = t278 * t284;
	t300 = t281 * t286;
	t299 = qJD(2) * t284;
	t298 = qJD(2) * t286;
	t296 = t278 * t298;
	t295 = t277 * t299;
	t294 = t280 * t298;
	t276 = sin(pkin(12));
	t279 = cos(pkin(12));
	t293 = -t279 * r_i_i_C(1) + t276 * r_i_i_C(2) - pkin(4);
	t292 = t276 * r_i_i_C(1) + t279 * r_i_i_C(2) + pkin(8) + qJ(4);
	t268 = -t277 * t301 + t280 * t286;
	t291 = t268 * t274 + t273 * t305;
	t290 = t281 * t273 + t274 * t302;
	t289 = t277 * t300 + t280 * t284;
	t285 = cos(qJ(3));
	t288 = -t285 * pkin(3) - t307 * t273 + t293 * t274 - pkin(2);
	t287 = t273 * qJD(5) + (-pkin(3) * t283 + t293 * t273 + t307 * t274) * qJD(3);
	t264 = -t281 * t295 + t294;
	t263 = t289 * qJD(2);
	t262 = t266 * qJD(2);
	t261 = -t281 * t294 + t295;
	t259 = t290 * qJD(3) + t273 * t296;
	t257 = t291 * qJD(3) - t263 * t273;
	t255 = -t308 * qJD(3) - t261 * t273;
	t1 = [0, t268 * qJD(4) - t292 * t263 + t288 * t264 - t287 * t289, t291 * qJD(5) + t307 * (-t263 * t274 + (-t268 * t273 + t274 * t305) * qJD(3)) + t293 * t257 + (t263 * t283 + (-t268 * t285 - t277 * t303) * qJD(3)) * pkin(3), t264, t257, 0; 0, t266 * qJD(4) - t292 * t261 + t288 * t262 + t287 * (-t277 * t284 + t280 * t300), -t308 * qJD(5) + t307 * (-t261 * t274 + (-t266 * t273 - t274 * t304) * qJD(3)) + t293 * t255 + (t261 * t283 + (-t266 * t285 + t280 * t303) * qJD(3)) * pkin(3), t262, t255, 0; 0, ((t288 * qJD(2) + qJD(4)) * t284 + (t292 * qJD(2) + t287) * t286) * t278, t290 * qJD(5) + t307 * (t274 * t296 + (-t273 * t302 + t274 * t281) * qJD(3)) + t293 * t259 + (-t283 * t296 + (-t281 * t283 - t285 * t302) * qJD(3)) * pkin(3), t278 * t299, t259, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:17
	% EndTime: 2019-10-09 22:07:18
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (508->83), mult. (905->142), div. (0->0), fcn. (904->14), ass. (0->54)
	t318 = qJ(3) + pkin(11);
	t314 = sin(t318);
	t316 = cos(t318);
	t324 = sin(qJ(3));
	t317 = pkin(12) + qJ(6);
	t313 = sin(t317);
	t315 = cos(t317);
	t339 = t313 * r_i_i_C(1) + t315 * r_i_i_C(2);
	t334 = qJD(6) * t339;
	t340 = t315 * r_i_i_C(1) - t313 * r_i_i_C(2);
	t337 = cos(pkin(12)) * pkin(5) + pkin(4) + t340;
	t356 = r_i_i_C(3) + pkin(9) + qJ(5);
	t328 = (t324 * pkin(3) + t337 * t314 - t356 * t316) * qJD(3) - t314 * qJD(5) + t316 * t334;
	t320 = sin(pkin(10));
	t325 = sin(qJ(2));
	t327 = cos(qJ(2));
	t354 = cos(pkin(10));
	t355 = cos(pkin(6));
	t338 = t355 * t354;
	t304 = t320 * t327 + t325 * t338;
	t321 = sin(pkin(6));
	t344 = t321 * t354;
	t294 = t304 * t316 - t314 * t344;
	t345 = t320 * t355;
	t306 = -t325 * t345 + t354 * t327;
	t352 = t320 * t321;
	t351 = t321 * t325;
	t350 = t321 * t327;
	t349 = qJD(2) * t325;
	t347 = qJD(2) * t350;
	t346 = t321 * t349;
	t336 = t327 * t338;
	t335 = -t306 * t314 + t316 * t352;
	t296 = t306 * t316 + t314 * t352;
	t333 = -t304 * t314 - t316 * t344;
	t332 = -t314 * t351 + t355 * t316;
	t298 = t355 * t314 + t316 * t351;
	t331 = sin(pkin(12)) * pkin(5) + qJ(4) + pkin(8) + t339;
	t330 = t340 * qJD(6) + qJD(4);
	t305 = t354 * t325 + t327 * t345;
	t326 = cos(qJ(3));
	t329 = -t326 * pkin(3) - t356 * t314 - t337 * t316 - pkin(2);
	t303 = t320 * t325 - t336;
	t302 = t306 * qJD(2);
	t301 = t305 * qJD(2);
	t300 = t304 * qJD(2);
	t299 = -qJD(2) * t336 + t320 * t349;
	t292 = t332 * qJD(3) + t316 * t347;
	t291 = t298 * qJD(3) + t314 * t347;
	t290 = t335 * qJD(3) - t301 * t316;
	t289 = t296 * qJD(3) - t301 * t314;
	t288 = t333 * qJD(3) - t299 * t316;
	t287 = t294 * qJD(3) - t299 * t314;
	t1 = [0, -t331 * t301 + t329 * t302 + t328 * t305 + t330 * t306, t296 * qJD(5) + t356 * t290 - t335 * t334 - t337 * t289 + (t301 * t324 + (-t306 * t326 - t324 * t352) * qJD(3)) * pkin(3), t302, t289, (-t290 * t313 + t302 * t315) * r_i_i_C(1) + (-t290 * t315 - t302 * t313) * r_i_i_C(2) + ((-t296 * t315 - t305 * t313) * r_i_i_C(1) + (t296 * t313 - t305 * t315) * r_i_i_C(2)) * qJD(6); 0, -t331 * t299 + t329 * t300 + t328 * t303 + t330 * t304, t294 * qJD(5) + t356 * t288 - t333 * t334 - t337 * t287 + (t299 * t324 + (-t304 * t326 + t324 * t344) * qJD(3)) * pkin(3), t300, t287, (-t288 * t313 + t300 * t315) * r_i_i_C(1) + (-t288 * t315 - t300 * t313) * r_i_i_C(2) + ((-t294 * t315 - t303 * t313) * r_i_i_C(1) + (t294 * t313 - t303 * t315) * r_i_i_C(2)) * qJD(6); 0, ((t329 * qJD(2) + t330) * t325 + (t331 * qJD(2) - t328) * t327) * t321, t298 * qJD(5) + t356 * t292 - t332 * t334 - t337 * t291 + (-t324 * t347 + (-t355 * t324 - t326 * t351) * qJD(3)) * pkin(3), t346, t291, (-t292 * t313 + t315 * t346) * r_i_i_C(1) + (-t292 * t315 - t313 * t346) * r_i_i_C(2) + ((-t298 * t315 + t313 * t350) * r_i_i_C(1) + (t298 * t313 + t315 * t350) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end