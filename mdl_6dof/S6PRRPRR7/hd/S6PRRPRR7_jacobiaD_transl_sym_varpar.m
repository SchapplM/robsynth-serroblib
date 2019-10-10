% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:19
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:37:20
	% EndTime: 2019-10-09 22:37:20
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
	t221 = sin(pkin(11));
	t223 = cos(pkin(11));
	t226 = sin(qJ(2));
	t224 = cos(pkin(6));
	t228 = cos(qJ(2));
	t238 = t224 * t228;
	t247 = -t221 * t226 + t223 * t238;
	t225 = sin(qJ(3));
	t227 = cos(qJ(3));
	t243 = r_i_i_C(3) + qJ(4);
	t245 = pkin(3) - r_i_i_C(2);
	t246 = t243 * t225 + t245 * t227 + pkin(2);
	t244 = pkin(8) + r_i_i_C(1);
	t222 = sin(pkin(6));
	t241 = t222 * t225;
	t240 = t222 * t227;
	t239 = t224 * t226;
	t236 = qJD(2) * t222 * t228;
	t217 = t221 * t228 + t223 * t239;
	t235 = -t217 * t227 + t223 * t241;
	t232 = t221 * t239 - t223 * t228;
	t234 = t221 * t241 - t227 * t232;
	t233 = t221 * t238 + t223 * t226;
	t231 = t224 * t225 + t226 * t240;
	t230 = qJD(2) * t246;
	t229 = qJD(4) * t225 + (-t245 * t225 + t243 * t227) * qJD(3);
	t214 = t233 * qJD(2);
	t212 = t247 * qJD(2);
	t210 = t231 * qJD(3) + t225 * t236;
	t208 = t234 * qJD(3) - t214 * t225;
	t206 = -t235 * qJD(3) + t212 * t225;
	t1 = [0, -t244 * t214 - t229 * t233 + t232 * t230, t234 * qJD(4) + t243 * (-t214 * t227 + (t221 * t240 + t225 * t232) * qJD(3)) - t245 * t208, t208, 0, 0; 0, t244 * t212 - t217 * t230 + t229 * t247, -t235 * qJD(4) + t243 * (t212 * t227 + (-t217 * t225 - t223 * t240) * qJD(3)) - t245 * t206, t206, 0, 0; 0, (t229 * t228 + (-t246 * t226 + t244 * t228) * qJD(2)) * t222, t231 * qJD(4) + t243 * (t227 * t236 + (t224 * t227 - t226 * t241) * qJD(3)) - t245 * t210, t210, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:20
	% EndTime: 2019-10-09 22:37:21
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (275->63), mult. (880->119), div. (0->0), fcn. (880->10), ass. (0->47)
	t302 = sin(qJ(3));
	t305 = cos(qJ(3));
	t301 = sin(qJ(5));
	t304 = cos(qJ(5));
	t320 = t304 * r_i_i_C(1) - t301 * r_i_i_C(2);
	t309 = t320 * qJD(5) + qJD(4);
	t319 = -r_i_i_C(1) * t301 - r_i_i_C(2) * t304;
	t317 = qJ(4) - t319;
	t328 = pkin(3) + pkin(9) + r_i_i_C(3);
	t307 = (t328 * t302 - t317 * t305) * qJD(3) - t309 * t302;
	t299 = sin(pkin(11));
	t303 = sin(qJ(2));
	t306 = cos(qJ(2));
	t334 = cos(pkin(11));
	t335 = cos(pkin(6));
	t318 = t335 * t334;
	t290 = t299 * t306 + t303 * t318;
	t300 = sin(pkin(6));
	t324 = t300 * t334;
	t337 = t290 * t305 - t302 * t324;
	t325 = t299 * t335;
	t292 = -t303 * t325 + t334 * t306;
	t332 = t300 * t302;
	t331 = t300 * t305;
	t330 = t300 * t306;
	t329 = qJD(2) * t303;
	t327 = t300 * t329;
	t326 = qJD(2) * t330;
	t316 = t306 * t318;
	t315 = pkin(4) + pkin(8) + t320;
	t314 = -t292 * t302 + t299 * t331;
	t313 = t292 * t305 + t299 * t332;
	t312 = t319 * qJD(5);
	t293 = t303 * t332 - t335 * t305;
	t311 = t335 * t302 + t303 * t331;
	t310 = -t290 * t302 - t305 * t324;
	t291 = t334 * t303 + t306 * t325;
	t308 = -t317 * t302 - t328 * t305 - pkin(2);
	t289 = t299 * t303 - t316;
	t288 = t292 * qJD(2);
	t287 = t291 * qJD(2);
	t286 = t290 * qJD(2);
	t285 = -qJD(2) * t316 + t299 * t329;
	t283 = t311 * qJD(3) + t302 * t326;
	t277 = t313 * qJD(3) - t287 * t302;
	t275 = t337 * qJD(3) - t285 * t302;
	t1 = [0, -t315 * t287 + t308 * t288 + t307 * t291 + t292 * t312, t309 * t313 + t317 * (t314 * qJD(3) - t287 * t305) - t328 * t277, t277, (t277 * t304 - t288 * t301) * r_i_i_C(1) + (-t277 * t301 - t288 * t304) * r_i_i_C(2) + ((-t291 * t304 + t301 * t314) * r_i_i_C(1) + (t291 * t301 + t304 * t314) * r_i_i_C(2)) * qJD(5), 0; 0, -t315 * t285 + t308 * t286 + t307 * t289 + t290 * t312, t309 * t337 + t317 * (t310 * qJD(3) - t285 * t305) - t328 * t275, t275, (t275 * t304 - t286 * t301) * r_i_i_C(1) + (-t275 * t301 - t286 * t304) * r_i_i_C(2) + ((-t289 * t304 + t301 * t310) * r_i_i_C(1) + (t289 * t301 + t304 * t310) * r_i_i_C(2)) * qJD(5), 0; 0, ((t308 * qJD(2) + t312) * t303 + (t315 * qJD(2) - t307) * t306) * t300, t309 * t311 + t317 * (-t293 * qJD(3) + t305 * t326) - t328 * t283, t283, (t283 * t304 - t301 * t327) * r_i_i_C(1) + (-t283 * t301 - t304 * t327) * r_i_i_C(2) + ((-t293 * t301 + t304 * t330) * r_i_i_C(1) + (-t293 * t304 - t301 * t330) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:20
	% EndTime: 2019-10-09 22:37:21
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (550->76), mult. (1306->127), div. (0->0), fcn. (1310->12), ass. (0->61)
	t343 = sin(qJ(3));
	t346 = cos(qJ(3));
	t338 = qJD(5) + qJD(6);
	t345 = cos(qJ(5));
	t339 = qJ(5) + qJ(6);
	t336 = sin(t339);
	t337 = cos(t339);
	t364 = t337 * r_i_i_C(1) - t336 * r_i_i_C(2);
	t388 = qJD(5) * pkin(5);
	t351 = t364 * t338 + t345 * t388 + qJD(4);
	t342 = sin(qJ(5));
	t363 = -t336 * r_i_i_C(1) - t337 * r_i_i_C(2);
	t353 = pkin(5) * t342 + qJ(4) - t363;
	t377 = pkin(3) + r_i_i_C(3) + pkin(10) + pkin(9);
	t349 = (t377 * t343 - t353 * t346) * qJD(3) - t351 * t343;
	t340 = sin(pkin(11));
	t344 = sin(qJ(2));
	t347 = cos(qJ(2));
	t386 = cos(pkin(11));
	t387 = cos(pkin(6));
	t361 = t387 * t386;
	t326 = t340 * t347 + t344 * t361;
	t341 = sin(pkin(6));
	t373 = t341 * t386;
	t390 = t326 * t346 - t343 * t373;
	t374 = t340 * t387;
	t328 = -t344 * t374 + t386 * t347;
	t384 = t341 * t343;
	t383 = t341 * t346;
	t382 = t341 * t347;
	t322 = t326 * qJD(2);
	t354 = -t326 * t343 - t346 * t373;
	t369 = t338 * t354 - t322;
	t360 = t347 * t361;
	t378 = qJD(2) * t344;
	t321 = -qJD(2) * t360 + t340 * t378;
	t311 = t390 * qJD(3) - t321 * t343;
	t325 = t340 * t344 - t360;
	t371 = t325 * t338 - t311;
	t381 = (t369 * t336 - t371 * t337) * r_i_i_C(1) + (t371 * t336 + t369 * t337) * r_i_i_C(2);
	t324 = t328 * qJD(2);
	t359 = -t328 * t343 + t340 * t383;
	t368 = t338 * t359 - t324;
	t327 = t386 * t344 + t347 * t374;
	t323 = t327 * qJD(2);
	t358 = t328 * t346 + t340 * t384;
	t313 = t358 * qJD(3) - t323 * t343;
	t370 = t327 * t338 - t313;
	t380 = (t368 * t336 - t370 * t337) * r_i_i_C(1) + (t370 * t336 + t368 * t337) * r_i_i_C(2);
	t329 = t344 * t384 - t387 * t346;
	t376 = t341 * t378;
	t357 = -t329 * t338 - t376;
	t355 = t387 * t343 + t344 * t383;
	t375 = qJD(2) * t382;
	t319 = t355 * qJD(3) + t343 * t375;
	t362 = t338 * t382 + t319;
	t379 = (t357 * t336 + t362 * t337) * r_i_i_C(1) + (-t362 * t336 + t357 * t337) * r_i_i_C(2);
	t356 = t345 * pkin(5) + pkin(4) + pkin(8) + t364;
	t352 = t363 * t338 - t342 * t388;
	t350 = -t353 * t343 - t377 * t346 - pkin(2);
	t1 = [0, -t356 * t323 + t350 * t324 + t349 * t327 + t352 * t328, -t377 * t313 + t351 * t358 + t353 * (t359 * qJD(3) - t323 * t346), t313, (t313 * t345 - t324 * t342 + (-t327 * t345 + t342 * t359) * qJD(5)) * pkin(5) + t380, t380; 0, -t356 * t321 + t350 * t322 + t349 * t325 + t352 * t326, -t377 * t311 + t351 * t390 + t353 * (t354 * qJD(3) - t321 * t346), t311, (t311 * t345 - t322 * t342 + (-t325 * t345 + t342 * t354) * qJD(5)) * pkin(5) + t381, t381; 0, ((t350 * qJD(2) + t352) * t344 + (t356 * qJD(2) - t349) * t347) * t341, -t377 * t319 + t351 * t355 + t353 * (-t329 * qJD(3) + t346 * t375), t319, (-t342 * t376 + t319 * t345 + (-t329 * t342 + t345 * t382) * qJD(5)) * pkin(5) + t379, t379;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end