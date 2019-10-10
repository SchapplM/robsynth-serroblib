% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 22:50:25
	% EndTime: 2019-10-09 22:50:25
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
	% StartTime: 2019-10-09 22:50:25
	% EndTime: 2019-10-09 22:50:25
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (180->39), mult. (333->75), div. (0->0), fcn. (304->10), ass. (0->38)
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t229 = sin(qJ(2));
	t227 = cos(pkin(6));
	t231 = cos(qJ(2));
	t247 = t227 * t231;
	t256 = -t224 * t229 + t226 * t247;
	t255 = r_i_i_C(3) + pkin(9) + pkin(8);
	t223 = qJ(3) + qJ(4);
	t220 = sin(t223);
	t222 = qJD(3) + qJD(4);
	t254 = t220 * t222;
	t221 = cos(t223);
	t253 = t221 * t222;
	t225 = sin(pkin(6));
	t252 = t222 * t225;
	t228 = sin(qJ(3));
	t250 = t225 * t228;
	t249 = t225 * t229;
	t248 = t227 * t229;
	t215 = t224 * t231 + t226 * t248;
	t210 = t256 * qJD(2);
	t240 = t226 * t252 - t210;
	t246 = (-t215 * t253 + t240 * t220) * r_i_i_C(1) + (t215 * t254 + t240 * t221) * r_i_i_C(2);
	t236 = t224 * t248 - t226 * t231;
	t237 = t224 * t247 + t226 * t229;
	t212 = t237 * qJD(2);
	t239 = -t224 * t252 + t212;
	t245 = (t239 * t220 + t236 * t253) * r_i_i_C(1) + (t239 * t221 - t236 * t254) * r_i_i_C(2);
	t241 = qJD(2) * t225 * t231;
	t235 = -t222 * t227 - t241;
	t243 = t222 * t249;
	t244 = (t235 * t220 - t221 * t243) * r_i_i_C(1) + (t220 * t243 + t235 * t221) * r_i_i_C(2);
	t230 = cos(qJ(3));
	t238 = t230 * pkin(3) + r_i_i_C(1) * t221 - r_i_i_C(2) * t220 + pkin(2);
	t234 = qJD(2) * t238;
	t233 = -pkin(3) * qJD(3) * t228 + (-r_i_i_C(1) * t220 - r_i_i_C(2) * t221) * t222;
	t1 = [0, -t255 * t212 - t233 * t237 + t236 * t234, (t212 * t228 + (-t224 * t250 + t230 * t236) * qJD(3)) * pkin(3) + t245, t245, 0, 0; 0, t255 * t210 - t215 * t234 + t233 * t256, (-t210 * t228 + (-t215 * t230 + t226 * t250) * qJD(3)) * pkin(3) + t246, t246, 0, 0; 0, (t233 * t231 + (-t238 * t229 + t255 * t231) * qJD(2)) * t225, (-t228 * t241 + (-t227 * t228 - t230 * t249) * qJD(3)) * pkin(3) + t244, t244, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:26
	% EndTime: 2019-10-09 22:50:26
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (391->49), mult. (642->86), div. (0->0), fcn. (615->10), ass. (0->45)
	t321 = pkin(4) - r_i_i_C(2);
	t319 = r_i_i_C(3) + qJ(5);
	t288 = sin(pkin(11));
	t290 = cos(pkin(11));
	t291 = cos(pkin(6));
	t295 = cos(qJ(2));
	t311 = t291 * t295;
	t307 = t290 * t311;
	t293 = sin(qJ(2));
	t310 = qJD(2) * t293;
	t272 = -qJD(2) * t307 + t288 * t310;
	t286 = qJD(3) + qJD(4);
	t289 = sin(pkin(6));
	t315 = t289 * t290;
	t323 = -t286 * t315 - t272;
	t287 = qJ(3) + qJ(4);
	t284 = sin(t287);
	t285 = cos(t287);
	t294 = cos(qJ(3));
	t322 = t294 * pkin(3) + t319 * t284 + t321 * t285 + pkin(2);
	t320 = r_i_i_C(1) + pkin(9) + pkin(8);
	t318 = t284 * t286;
	t317 = t285 * t286;
	t316 = t288 * t289;
	t292 = sin(qJ(3));
	t314 = t289 * t292;
	t313 = t289 * t293;
	t312 = t291 * t293;
	t308 = t285 * t313;
	t306 = qJD(2) * t289 * t295;
	t304 = t288 * t311 + t290 * t293;
	t274 = t304 * qJD(2);
	t305 = t286 * t316 - t274;
	t277 = t288 * t295 + t290 * t312;
	t303 = t288 * t312 - t290 * t295;
	t302 = t286 * t291 + t306;
	t258 = t277 * t317 + t323 * t284;
	t301 = -(-t277 * t285 + t284 * t315) * qJD(5) + t319 * (-t277 * t318 + t323 * t285) - t321 * t258;
	t260 = t305 * t284 - t303 * t317;
	t300 = -(-t284 * t316 + t285 * t303) * qJD(5) + t319 * (t305 * t285 + t303 * t318) - t321 * t260;
	t265 = t302 * t284 + t286 * t308;
	t299 = -(-t291 * t284 - t308) * qJD(5) + t319 * (t302 * t285 - t313 * t318) - t321 * t265;
	t298 = qJD(2) * t322;
	t297 = -pkin(3) * qJD(3) * t292 + qJD(5) * t284 + (-t321 * t284 + t319 * t285) * t286;
	t1 = [0, -t320 * t274 - t297 * t304 + t303 * t298, (t274 * t292 + (-t288 * t314 + t294 * t303) * qJD(3)) * pkin(3) + t300, t300, t260, 0; 0, -t320 * t272 - t277 * t298 + t297 * (-t288 * t293 + t307), (t272 * t292 + (-t277 * t294 + t290 * t314) * qJD(3)) * pkin(3) + t301, t301, t258, 0; 0, (-t322 * t310 + (t320 * qJD(2) + t297) * t295) * t289, (-t292 * t306 + (-t291 * t292 - t294 * t313) * qJD(3)) * pkin(3) + t299, t299, t265, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:27
	% EndTime: 2019-10-09 22:50:27
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (734->80), mult. (1250->142), div. (0->0), fcn. (1245->12), ass. (0->59)
	t383 = sin(qJ(6));
	t386 = cos(qJ(6));
	t404 = t386 * r_i_i_C(1) - t383 * r_i_i_C(2);
	t433 = t404 * qJD(6) + qJD(5);
	t403 = -t383 * r_i_i_C(1) - t386 * r_i_i_C(2);
	t434 = qJ(5) - t403;
	t417 = pkin(4) + pkin(10) + r_i_i_C(3);
	t380 = qJ(3) + qJ(4);
	t377 = sin(t380);
	t379 = qJD(3) + qJD(4);
	t384 = sin(qJ(3));
	t378 = cos(t380);
	t422 = t378 * t379;
	t390 = -(-t417 * t379 + t433) * t377 + qJD(3) * t384 * pkin(3) - t434 * t422;
	t385 = sin(qJ(2));
	t388 = cos(qJ(2));
	t381 = sin(pkin(11));
	t425 = cos(pkin(6));
	t410 = t381 * t425;
	t424 = cos(pkin(11));
	t368 = -t385 * t410 + t424 * t388;
	t423 = t377 * t379;
	t382 = sin(pkin(6));
	t421 = t381 * t382;
	t420 = t382 * t385;
	t419 = t382 * t388;
	t418 = qJD(2) * t385;
	t414 = t377 * t420;
	t413 = t378 * t420;
	t412 = t382 * t418;
	t411 = qJD(2) * t419;
	t409 = t382 * t424;
	t406 = t377 * t409;
	t405 = t378 * t409;
	t367 = t424 * t385 + t388 * t410;
	t363 = t367 * qJD(2);
	t401 = t379 * t421 - t363;
	t400 = t425 * t424;
	t398 = t388 * t400;
	t397 = t403 * qJD(6);
	t396 = pkin(5) + pkin(9) + pkin(8) + t404;
	t395 = t425 * t379 + t411;
	t366 = t381 * t388 + t385 * t400;
	t387 = cos(qJ(3));
	t394 = -t387 * pkin(3) - t377 * t434 - t417 * t378 - pkin(2);
	t361 = -qJD(2) * t398 + t381 * t418;
	t340 = -t361 * t377 + t366 * t422 - t379 * t406;
	t393 = t433 * (t366 * t378 - t406) + t434 * (-t361 * t378 - t366 * t423 - t379 * t405) - t417 * t340;
	t342 = t368 * t422 + t401 * t377;
	t392 = t433 * (t368 * t378 + t377 * t421) + t434 * (-t368 * t423 + t401 * t378) - t417 * t342;
	t348 = t395 * t377 + t379 * t413;
	t391 = t433 * (t425 * t377 + t413) + t434 * (t395 * t378 - t379 * t414) - t417 * t348;
	t365 = t381 * t385 - t398;
	t364 = t368 * qJD(2);
	t362 = t366 * qJD(2);
	t359 = -t425 * t378 + t414;
	t354 = t368 * t377 - t378 * t421;
	t352 = t366 * t377 + t405;
	t1 = [0, -t396 * t363 + t394 * t364 + t390 * t367 + t368 * t397, (t363 * t384 + (-t368 * t387 - t384 * t421) * qJD(3)) * pkin(3) + t392, t392, t342, (t342 * t386 - t364 * t383) * r_i_i_C(1) + (-t342 * t383 - t364 * t386) * r_i_i_C(2) + ((-t354 * t383 - t367 * t386) * r_i_i_C(1) + (-t354 * t386 + t367 * t383) * r_i_i_C(2)) * qJD(6); 0, -t396 * t361 + t394 * t362 + t390 * t365 + t366 * t397, (t361 * t384 + (-t366 * t387 + t384 * t409) * qJD(3)) * pkin(3) + t393, t393, t340, (t340 * t386 - t362 * t383) * r_i_i_C(1) + (-t340 * t383 - t362 * t386) * r_i_i_C(2) + ((-t352 * t383 - t365 * t386) * r_i_i_C(1) + (-t352 * t386 + t365 * t383) * r_i_i_C(2)) * qJD(6); 0, ((t394 * qJD(2) + t397) * t385 + (t396 * qJD(2) - t390) * t388) * t382, (-t384 * t411 + (-t425 * t384 - t387 * t420) * qJD(3)) * pkin(3) + t391, t391, t348, (t348 * t386 - t383 * t412) * r_i_i_C(1) + (-t348 * t383 - t386 * t412) * r_i_i_C(2) + ((-t359 * t383 + t386 * t419) * r_i_i_C(1) + (-t359 * t386 - t383 * t419) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end