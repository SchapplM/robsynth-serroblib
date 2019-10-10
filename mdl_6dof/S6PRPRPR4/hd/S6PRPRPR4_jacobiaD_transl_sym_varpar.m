% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
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
	% StartTime: 2019-10-09 21:35:32
	% EndTime: 2019-10-09 21:35:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->14), mult. (103->28), div. (0->0), fcn. (90->8), ass. (0->16)
	t153 = sin(pkin(10));
	t156 = cos(pkin(10));
	t159 = cos(qJ(2));
	t157 = cos(pkin(6));
	t158 = sin(qJ(2));
	t165 = t157 * t158;
	t168 = t153 * t165 - t156 * t159;
	t167 = r_i_i_C(3) + qJ(3);
	t164 = t157 * t159;
	t162 = qJD(2) * t167;
	t161 = -cos(pkin(11)) * r_i_i_C(1) + sin(pkin(11)) * r_i_i_C(2) - pkin(2);
	t160 = t153 * t159 + t156 * t165;
	t154 = sin(pkin(6));
	t150 = t168 * qJD(2);
	t148 = t160 * qJD(2);
	t1 = [0, -t168 * qJD(3) - (t153 * t164 + t156 * t158) * t162 - t161 * t150, -t150, 0, 0, 0; 0, t160 * qJD(3) - (t153 * t158 - t156 * t164) * t162 + t161 * t148, t148, 0, 0, 0; 0, (qJD(3) * t158 + (t161 * t158 + t167 * t159) * qJD(2)) * t154, t154 * qJD(2) * t158, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:32
	% EndTime: 2019-10-09 21:35:32
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (98->33), mult. (216->67), div. (0->0), fcn. (199->9), ass. (0->30)
	t212 = r_i_i_C(3) + pkin(8) + qJ(3);
	t192 = sin(pkin(10));
	t193 = sin(pkin(6));
	t211 = t192 * t193;
	t194 = cos(pkin(10));
	t210 = t193 * t194;
	t197 = sin(qJ(2));
	t209 = t193 * t197;
	t195 = cos(pkin(6));
	t208 = t195 * t197;
	t198 = cos(qJ(2));
	t207 = t195 * t198;
	t206 = qJD(2) * t197;
	t205 = qJD(2) * t198;
	t204 = t192 * t206;
	t203 = t194 * t205;
	t191 = pkin(11) + qJ(4);
	t189 = sin(t191);
	t190 = cos(t191);
	t202 = t189 * r_i_i_C(1) + t190 * r_i_i_C(2);
	t201 = -t190 * r_i_i_C(1) + t189 * r_i_i_C(2) - cos(pkin(11)) * pkin(3) - pkin(2);
	t183 = t192 * t198 + t194 * t208;
	t200 = t192 * t207 + t194 * t197;
	t199 = qJD(4) * t202;
	t185 = -t192 * t208 + t194 * t198;
	t181 = -t195 * t204 + t203;
	t180 = t200 * qJD(2);
	t179 = t183 * qJD(2);
	t178 = -t195 * t203 + t204;
	t1 = [0, t185 * qJD(3) - t212 * t180 + t201 * t181 + t200 * t199, t181, t202 * t180 + ((-t185 * t190 - t189 * t211) * r_i_i_C(1) + (t185 * t189 - t190 * t211) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t183 * qJD(3) - t212 * t178 - (-t192 * t197 + t194 * t207) * t199 + t201 * t179, t179, t202 * t178 + ((-t183 * t190 + t189 * t210) * r_i_i_C(1) + (t183 * t189 + t190 * t210) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (qJD(3) * t197 - t198 * t199 + (t201 * t197 + t212 * t198) * qJD(2)) * t193, t193 * t206, -t202 * t193 * t205 + ((-t189 * t195 - t190 * t209) * r_i_i_C(1) + (t189 * t209 - t190 * t195) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:33
	% EndTime: 2019-10-09 21:35:33
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (258->45), mult. (529->81), div. (0->0), fcn. (515->11), ass. (0->40)
	t274 = sin(pkin(10));
	t277 = cos(pkin(10));
	t281 = cos(qJ(2));
	t278 = cos(pkin(6));
	t280 = sin(qJ(2));
	t296 = t278 * t280;
	t263 = t274 * t281 + t277 * t296;
	t272 = pkin(11) + qJ(4);
	t270 = sin(t272);
	t271 = cos(t272);
	t275 = sin(pkin(6));
	t298 = t275 * t277;
	t302 = -t263 * t271 + t270 * t298;
	t301 = r_i_i_C(3) + qJ(5);
	t299 = t274 * t275;
	t297 = t275 * t280;
	t295 = t278 * t281;
	t294 = qJD(2) * t280;
	t293 = qJD(2) * t281;
	t291 = t275 * t293;
	t290 = t274 * t294;
	t289 = t277 * t293;
	t273 = sin(pkin(12));
	t276 = cos(pkin(12));
	t288 = -t276 * r_i_i_C(1) + t273 * r_i_i_C(2) - pkin(4);
	t287 = t273 * r_i_i_C(1) + t276 * r_i_i_C(2) + pkin(8) + qJ(3);
	t265 = -t274 * t296 + t277 * t281;
	t286 = t265 * t271 + t270 * t299;
	t285 = t278 * t270 + t271 * t297;
	t284 = t274 * t295 + t277 * t280;
	t283 = -t301 * t270 + t288 * t271 - cos(pkin(11)) * pkin(3) - pkin(2);
	t282 = t270 * qJD(5) + (t288 * t270 + t301 * t271) * qJD(4);
	t261 = -t278 * t290 + t289;
	t260 = t284 * qJD(2);
	t259 = t263 * qJD(2);
	t258 = -t278 * t289 + t290;
	t256 = t285 * qJD(4) + t270 * t291;
	t254 = t286 * qJD(4) - t260 * t270;
	t252 = -t302 * qJD(4) - t258 * t270;
	t1 = [0, t265 * qJD(3) - t287 * t260 + t283 * t261 - t282 * t284, t261, t286 * qJD(5) + t301 * (-t260 * t271 + (-t265 * t270 + t271 * t299) * qJD(4)) + t288 * t254, t254, 0; 0, t263 * qJD(3) - t287 * t258 + t282 * (-t274 * t280 + t277 * t295) + t283 * t259, t259, -t302 * qJD(5) + t301 * (-t258 * t271 + (-t263 * t270 - t271 * t298) * qJD(4)) + t288 * t252, t252, 0; 0, (t280 * qJD(3) + t282 * t281 + (t283 * t280 + t287 * t281) * qJD(2)) * t275, t275 * t294, t285 * qJD(5) + t301 * (t271 * t291 + (-t270 * t297 + t271 * t278) * qJD(4)) + t288 * t256, t256, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:33
	% EndTime: 2019-10-09 21:35:33
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (490->73), mult. (844->126), div. (0->0), fcn. (851->13), ass. (0->52)
	t315 = pkin(11) + qJ(4);
	t311 = sin(t315);
	t313 = cos(t315);
	t314 = pkin(12) + qJ(6);
	t310 = sin(t314);
	t312 = cos(t314);
	t334 = t310 * r_i_i_C(1) + t312 * r_i_i_C(2);
	t329 = qJD(6) * t334;
	t335 = t312 * r_i_i_C(1) - t310 * r_i_i_C(2);
	t332 = cos(pkin(12)) * pkin(5) + pkin(4) + t335;
	t351 = r_i_i_C(3) + pkin(9) + qJ(5);
	t323 = (t332 * t311 - t351 * t313) * qJD(4) - t311 * qJD(5) + t313 * t329;
	t317 = sin(pkin(10));
	t321 = sin(qJ(2));
	t322 = cos(qJ(2));
	t349 = cos(pkin(10));
	t350 = cos(pkin(6));
	t333 = t350 * t349;
	t301 = t317 * t322 + t321 * t333;
	t318 = sin(pkin(6));
	t339 = t318 * t349;
	t291 = t301 * t313 - t311 * t339;
	t340 = t317 * t350;
	t303 = -t321 * t340 + t349 * t322;
	t347 = t317 * t318;
	t346 = t318 * t321;
	t345 = t318 * t322;
	t344 = qJD(2) * t321;
	t342 = qJD(2) * t345;
	t341 = t318 * t344;
	t331 = t322 * t333;
	t330 = -t303 * t311 + t313 * t347;
	t293 = t303 * t313 + t311 * t347;
	t328 = -t301 * t311 - t313 * t339;
	t327 = -t311 * t346 + t350 * t313;
	t295 = t350 * t311 + t313 * t346;
	t326 = sin(pkin(12)) * pkin(5) + pkin(8) + qJ(3) + t334;
	t325 = t335 * qJD(6) + qJD(3);
	t302 = t349 * t321 + t322 * t340;
	t324 = -t351 * t311 - t332 * t313 - cos(pkin(11)) * pkin(3) - pkin(2);
	t300 = t317 * t321 - t331;
	t299 = t303 * qJD(2);
	t298 = t302 * qJD(2);
	t297 = t301 * qJD(2);
	t296 = -qJD(2) * t331 + t317 * t344;
	t289 = t327 * qJD(4) + t313 * t342;
	t288 = t295 * qJD(4) + t311 * t342;
	t287 = t330 * qJD(4) - t298 * t313;
	t286 = t293 * qJD(4) - t298 * t311;
	t285 = t328 * qJD(4) - t296 * t313;
	t284 = t291 * qJD(4) - t296 * t311;
	t1 = [0, -t326 * t298 + t324 * t299 + t323 * t302 + t325 * t303, t299, t293 * qJD(5) - t332 * t286 + t351 * t287 - t330 * t329, t286, (-t287 * t310 + t299 * t312) * r_i_i_C(1) + (-t287 * t312 - t299 * t310) * r_i_i_C(2) + ((-t293 * t312 - t302 * t310) * r_i_i_C(1) + (t293 * t310 - t302 * t312) * r_i_i_C(2)) * qJD(6); 0, -t326 * t296 + t324 * t297 + t323 * t300 + t325 * t301, t297, t291 * qJD(5) - t332 * t284 + t351 * t285 - t328 * t329, t284, (-t285 * t310 + t297 * t312) * r_i_i_C(1) + (-t285 * t312 - t297 * t310) * r_i_i_C(2) + ((-t291 * t312 - t300 * t310) * r_i_i_C(1) + (t291 * t310 - t300 * t312) * r_i_i_C(2)) * qJD(6); 0, ((t324 * qJD(2) + t325) * t321 + (t326 * qJD(2) - t323) * t322) * t318, t341, t295 * qJD(5) - t332 * t288 + t351 * t289 - t327 * t329, t288, (-t289 * t310 + t312 * t341) * r_i_i_C(1) + (-t289 * t312 - t310 * t341) * r_i_i_C(2) + ((-t295 * t312 + t310 * t345) * r_i_i_C(1) + (t295 * t310 + t312 * t345) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end