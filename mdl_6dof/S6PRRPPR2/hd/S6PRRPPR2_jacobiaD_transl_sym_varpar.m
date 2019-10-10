% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
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
	% StartTime: 2019-10-09 22:09:07
	% EndTime: 2019-10-09 22:09:07
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
	% StartTime: 2019-10-09 22:09:07
	% EndTime: 2019-10-09 22:09:07
	% DurationCPUTime: 0.20s
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
	% StartTime: 2019-10-09 22:09:07
	% EndTime: 2019-10-09 22:09:08
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (233->52), mult. (490->93), div. (0->0), fcn. (464->10), ass. (0->41)
	t267 = -pkin(4) + r_i_i_C(2);
	t266 = r_i_i_C(1) + qJ(4) + pkin(8);
	t265 = r_i_i_C(3) + qJ(5);
	t239 = sin(pkin(10));
	t240 = sin(pkin(6));
	t264 = t239 * t240;
	t241 = cos(pkin(10));
	t263 = t240 * t241;
	t244 = sin(qJ(3));
	t262 = t240 * t244;
	t245 = sin(qJ(2));
	t261 = t240 * t245;
	t242 = cos(pkin(6));
	t260 = t242 * t245;
	t247 = cos(qJ(2));
	t259 = t242 * t247;
	t258 = qJD(2) * t245;
	t257 = qJD(2) * t247;
	t256 = t240 * t257;
	t255 = t239 * t258;
	t254 = t241 * t257;
	t230 = t239 * t247 + t241 * t260;
	t238 = qJ(3) + pkin(11);
	t236 = sin(t238);
	t237 = cos(t238);
	t253 = -t230 * t237 + t236 * t263;
	t232 = -t239 * t260 + t241 * t247;
	t252 = t232 * t237 + t236 * t264;
	t251 = t242 * t236 + t237 * t261;
	t250 = t239 * t259 + t241 * t245;
	t246 = cos(qJ(3));
	t249 = -t246 * pkin(3) - t265 * t236 + t267 * t237 - pkin(2);
	t248 = t236 * qJD(5) + (-pkin(3) * t244 + t267 * t236 + t265 * t237) * qJD(3);
	t228 = -t242 * t255 + t254;
	t227 = t250 * qJD(2);
	t226 = t230 * qJD(2);
	t225 = -t242 * t254 + t255;
	t223 = t251 * qJD(3) + t236 * t256;
	t221 = t252 * qJD(3) - t227 * t236;
	t219 = -t253 * qJD(3) - t225 * t236;
	t1 = [0, t232 * qJD(4) - t266 * t227 + t249 * t228 - t248 * t250, t252 * qJD(5) + t265 * (-t227 * t237 + (-t232 * t236 + t237 * t264) * qJD(3)) + t267 * t221 + (t227 * t244 + (-t232 * t246 - t239 * t262) * qJD(3)) * pkin(3), t228, t221, 0; 0, t230 * qJD(4) - t266 * t225 + t249 * t226 + t248 * (-t239 * t245 + t241 * t259), -t253 * qJD(5) + t265 * (-t225 * t237 + (-t230 * t236 - t237 * t263) * qJD(3)) + t267 * t219 + (t225 * t244 + (-t230 * t246 + t241 * t262) * qJD(3)) * pkin(3), t226, t219, 0; 0, ((t249 * qJD(2) + qJD(4)) * t245 + (t266 * qJD(2) + t248) * t247) * t240, t251 * qJD(5) + t265 * (t237 * t256 + (-t236 * t261 + t237 * t242) * qJD(3)) + t267 * t223 + (-t244 * t256 + (-t242 * t244 - t246 * t261) * qJD(3)) * pkin(3), t240 * t258, t223, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:08
	% EndTime: 2019-10-09 22:09:08
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (462->76), mult. (966->136), div. (0->0), fcn. (960->12), ass. (0->50)
	t316 = qJ(3) + pkin(11);
	t314 = sin(t316);
	t315 = cos(t316);
	t321 = sin(qJ(3));
	t320 = sin(qJ(6));
	t323 = cos(qJ(6));
	t339 = t323 * r_i_i_C(1) - t320 * r_i_i_C(2);
	t328 = t339 * qJD(6) + qJD(5);
	t338 = -t320 * r_i_i_C(1) - t323 * r_i_i_C(2);
	t336 = qJ(5) - t338;
	t347 = pkin(4) + pkin(9) + r_i_i_C(3);
	t326 = (t321 * pkin(3) + t347 * t314 - t336 * t315) * qJD(3) - t328 * t314;
	t317 = sin(pkin(10));
	t322 = sin(qJ(2));
	t325 = cos(qJ(2));
	t353 = cos(pkin(10));
	t354 = cos(pkin(6));
	t337 = t354 * t353;
	t306 = t317 * t325 + t322 * t337;
	t318 = sin(pkin(6));
	t343 = t318 * t353;
	t357 = t306 * t315 - t314 * t343;
	t344 = t317 * t354;
	t308 = -t322 * t344 + t353 * t325;
	t351 = t317 * t318;
	t350 = t318 * t322;
	t349 = t318 * t325;
	t348 = qJD(2) * t322;
	t346 = qJD(2) * t349;
	t345 = t318 * t348;
	t335 = t325 * t337;
	t334 = -t308 * t314 + t315 * t351;
	t333 = t308 * t315 + t314 * t351;
	t332 = pkin(5) + qJ(4) + pkin(8) + t339;
	t299 = t314 * t350 - t354 * t315;
	t331 = t354 * t314 + t315 * t350;
	t330 = -t306 * t314 - t315 * t343;
	t329 = t338 * qJD(6) + qJD(4);
	t307 = t353 * t322 + t325 * t344;
	t324 = cos(qJ(3));
	t327 = -t324 * pkin(3) - t336 * t314 - t347 * t315 - pkin(2);
	t305 = t317 * t322 - t335;
	t304 = t308 * qJD(2);
	t303 = t307 * qJD(2);
	t302 = t306 * qJD(2);
	t301 = -qJD(2) * t335 + t317 * t348;
	t293 = t331 * qJD(3) + t314 * t346;
	t291 = t333 * qJD(3) - t303 * t314;
	t289 = t357 * qJD(3) - t301 * t314;
	t1 = [0, -t332 * t303 + t327 * t304 + t326 * t307 + t329 * t308, t328 * t333 + t336 * (t334 * qJD(3) - t303 * t315) - t347 * t291 + (t303 * t321 + (-t308 * t324 - t321 * t351) * qJD(3)) * pkin(3), t304, t291, (t291 * t323 - t304 * t320) * r_i_i_C(1) + (-t291 * t320 - t304 * t323) * r_i_i_C(2) + ((-t307 * t323 + t320 * t334) * r_i_i_C(1) + (t307 * t320 + t323 * t334) * r_i_i_C(2)) * qJD(6); 0, -t332 * t301 + t327 * t302 + t326 * t305 + t329 * t306, t328 * t357 + t336 * (t330 * qJD(3) - t301 * t315) - t347 * t289 + (t301 * t321 + (-t306 * t324 + t321 * t343) * qJD(3)) * pkin(3), t302, t289, (t289 * t323 - t302 * t320) * r_i_i_C(1) + (-t289 * t320 - t302 * t323) * r_i_i_C(2) + ((-t305 * t323 + t320 * t330) * r_i_i_C(1) + (t305 * t320 + t323 * t330) * r_i_i_C(2)) * qJD(6); 0, ((t327 * qJD(2) + t329) * t322 + (t332 * qJD(2) - t326) * t325) * t318, t328 * t331 + t336 * (-t299 * qJD(3) + t315 * t346) - t347 * t293 + (-t321 * t346 + (-t354 * t321 - t324 * t350) * qJD(3)) * pkin(3), t345, t293, (t293 * t323 - t320 * t345) * r_i_i_C(1) + (-t293 * t320 - t323 * t345) * r_i_i_C(2) + ((-t299 * t320 + t323 * t349) * r_i_i_C(1) + (-t299 * t323 - t320 * t349) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end