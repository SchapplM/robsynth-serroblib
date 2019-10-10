% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
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
	% StartTime: 2019-10-09 21:28:09
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
	t132 = sin(pkin(10));
	t134 = cos(pkin(10));
	t137 = cos(qJ(2));
	t135 = cos(pkin(6));
	t136 = sin(qJ(2));
	t142 = t135 * t136;
	t146 = t132 * t142 - t134 * t137;
	t145 = -pkin(2) - r_i_i_C(1);
	t144 = r_i_i_C(3) + qJ(3);
	t141 = t135 * t137;
	t139 = qJD(2) * t144;
	t138 = t132 * t137 + t134 * t142;
	t133 = sin(pkin(6));
	t130 = t146 * qJD(2);
	t128 = t138 * qJD(2);
	t1 = [0, -t146 * qJD(3) - t145 * t130 - (t132 * t141 + t134 * t136) * t139, -t130, 0, 0, 0; 0, t138 * qJD(3) + t145 * t128 - (t132 * t136 - t134 * t141) * t139, t128, 0, 0, 0; 0, (qJD(3) * t136 + (t145 * t136 + t144 * t137) * qJD(2)) * t133, t133 * qJD(2) * t136, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (39->16), mult. (133->30), div. (0->0), fcn. (120->8), ass. (0->20)
	t103 = cos(qJ(2));
	t97 = sin(pkin(10));
	t112 = t97 * t103;
	t100 = cos(pkin(10));
	t111 = t100 * t103;
	t101 = cos(pkin(6));
	t102 = sin(qJ(2));
	t110 = t101 * t102;
	t109 = qJD(2) * t102;
	t108 = t97 * t109;
	t107 = qJD(2) * t111;
	t96 = sin(pkin(11));
	t99 = cos(pkin(11));
	t106 = t96 * r_i_i_C(1) + t99 * r_i_i_C(2) + qJ(3);
	t105 = -t99 * r_i_i_C(1) + t96 * r_i_i_C(2) - pkin(2) - pkin(3);
	t104 = t100 * t110 + t112;
	t98 = sin(pkin(6));
	t93 = -t101 * t108 + t107;
	t91 = t104 * qJD(2);
	t1 = [0, -(t97 * t110 - t111) * qJD(3) - t106 * (t100 * t102 + t101 * t112) * qJD(2) + t105 * t93, t93, 0, 0, 0; 0, t104 * qJD(3) - t106 * (-t101 * t107 + t108) + t105 * t91, t91, 0, 0, 0; 0, (t102 * qJD(3) + (t105 * t102 + t106 * t103) * qJD(2)) * t98, t98 * t109, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:09
	% EndTime: 2019-10-09 21:28:10
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (132->49), mult. (436->94), div. (0->0), fcn. (430->10), ass. (0->40)
	t277 = -pkin(2) - pkin(3);
	t276 = pkin(8) + r_i_i_C(3);
	t253 = sin(pkin(6));
	t257 = sin(qJ(5));
	t275 = t253 * t257;
	t259 = cos(qJ(5));
	t274 = t253 * t259;
	t256 = cos(pkin(6));
	t258 = sin(qJ(2));
	t273 = t256 * t258;
	t260 = cos(qJ(2));
	t272 = t256 * t260;
	t271 = qJD(2) * t258;
	t270 = qJD(2) * t260;
	t252 = sin(pkin(10));
	t269 = t252 * t271;
	t268 = t253 * t271;
	t255 = cos(pkin(10));
	t267 = t255 * t270;
	t239 = -t256 * t267 + t269;
	t244 = t252 * t260 + t255 * t273;
	t240 = t244 * qJD(2);
	t251 = sin(pkin(11));
	t254 = cos(pkin(11));
	t266 = -t239 * t254 + t240 * t251;
	t245 = t252 * t272 + t255 * t258;
	t241 = t245 * qJD(2);
	t242 = -t256 * t269 + t267;
	t265 = -t241 * t254 + t242 * t251;
	t264 = t257 * r_i_i_C(1) + t259 * r_i_i_C(2);
	t263 = t251 * t260 - t254 * t258;
	t262 = t259 * r_i_i_C(1) - t257 * r_i_i_C(2) + pkin(4);
	t261 = qJD(5) * t264;
	t246 = -t252 * t273 + t255 * t260;
	t243 = t252 * t258 - t255 * t272;
	t238 = t263 * t253;
	t235 = -t253 * t254 * t270 - t251 * t268;
	t232 = t245 * t251 + t246 * t254;
	t230 = t243 * t251 + t244 * t254;
	t1 = [0, -t241 * qJ(3) + t246 * qJD(3) + t277 * t242 - t276 * t265 - (-t245 * t254 + t246 * t251) * t261 + t262 * (-t241 * t251 - t242 * t254), t242, 0, -t264 * t265 + ((-t232 * t259 + t252 * t275) * r_i_i_C(1) + (t232 * t257 + t252 * t274) * r_i_i_C(2)) * qJD(5), 0; 0, -t239 * qJ(3) + t244 * qJD(3) + t277 * t240 - t276 * t266 - (-t243 * t254 + t244 * t251) * t261 + t262 * (-t239 * t251 - t240 * t254), t240, 0, -t264 * t266 + ((-t230 * t259 - t255 * t275) * r_i_i_C(1) + (t230 * t257 - t255 * t274) * r_i_i_C(2)) * qJD(5), 0; 0, t276 * t235 + (-(t251 * t258 + t254 * t260) * t261 + qJD(3) * t258 + (qJ(3) * t260 + t277 * t258 + t262 * t263) * qJD(2)) * t253, t268, 0, t264 * t235 + ((t238 * t259 + t256 * t257) * r_i_i_C(1) + (-t238 * t257 + t256 * t259) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:11
	% EndTime: 2019-10-09 21:28:11
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (446->92), mult. (1402->170), div. (0->0), fcn. (1504->12), ass. (0->59)
	t393 = sin(pkin(11));
	t395 = sin(pkin(6));
	t396 = cos(pkin(11));
	t401 = sin(qJ(2));
	t404 = cos(qJ(2));
	t378 = (-t393 * t404 + t396 * t401) * t395;
	t427 = -pkin(3) - pkin(2);
	t426 = pkin(9) + r_i_i_C(3);
	t400 = sin(qJ(5));
	t424 = t395 * t400;
	t403 = cos(qJ(5));
	t423 = t395 * t403;
	t398 = cos(pkin(6));
	t421 = t398 * t401;
	t420 = t398 * t404;
	t419 = qJD(2) * t401;
	t418 = qJD(2) * t404;
	t399 = sin(qJ(6));
	t417 = qJD(6) * t399;
	t402 = cos(qJ(6));
	t416 = qJD(6) * t402;
	t394 = sin(pkin(10));
	t415 = t394 * t419;
	t414 = t395 * t419;
	t397 = cos(pkin(10));
	t413 = t397 * t418;
	t379 = -t398 * t413 + t415;
	t384 = t394 * t404 + t397 * t421;
	t380 = t384 * qJD(2);
	t351 = -t379 * t393 - t380 * t396;
	t352 = -t379 * t396 + t380 * t393;
	t385 = t394 * t420 + t397 * t401;
	t381 = t385 * qJD(2);
	t382 = -t398 * t415 + t413;
	t355 = -t381 * t393 - t382 * t396;
	t356 = -t381 * t396 + t382 * t393;
	t383 = t394 * t401 - t397 * t420;
	t361 = -t383 * t396 + t384 * t393;
	t386 = -t394 * t421 + t397 * t404;
	t365 = -t385 * t396 + t386 * t393;
	t368 = t378 * t403 - t398 * t400;
	t412 = -t378 * t400 - t398 * t403;
	t362 = t383 * t393 + t384 * t396;
	t366 = t385 * t393 + t386 * t396;
	t411 = t402 * r_i_i_C(1) - t399 * r_i_i_C(2) + pkin(5);
	t410 = -t362 * t400 + t397 * t423;
	t348 = t362 * t403 + t397 * t424;
	t409 = -t366 * t400 - t394 * t423;
	t408 = -t366 * t403 + t394 * t424;
	t407 = qJD(6) * (-t399 * r_i_i_C(1) - t402 * r_i_i_C(2));
	t406 = t426 * t400 + t411 * t403 + pkin(4);
	t405 = t403 * t407 + (-t411 * t400 + t426 * t403) * qJD(5);
	t377 = (t393 * t401 + t396 * t404) * t395;
	t374 = qJD(2) * t378;
	t373 = -t395 * t396 * t418 - t393 * t414;
	t346 = t412 * qJD(5) - t373 * t403;
	t344 = t409 * qJD(5) + t356 * t403;
	t342 = t410 * qJD(5) + t352 * t403;
	t1 = [0, (-t356 * t399 - t366 * t416) * r_i_i_C(1) + (-t356 * t402 + t366 * t417) * r_i_i_C(2) - t356 * pkin(8) - t381 * qJ(3) + t386 * qJD(3) + t427 * t382 + t406 * t355 + t405 * t365, t382, 0, t426 * t344 + t409 * t407 + t411 * (t408 * qJD(5) - t356 * t400), (-t344 * t399 + t355 * t402) * r_i_i_C(1) + (-t344 * t402 - t355 * t399) * r_i_i_C(2) + ((-t365 * t399 + t402 * t408) * r_i_i_C(1) + (-t365 * t402 - t399 * t408) * r_i_i_C(2)) * qJD(6); 0, (-t352 * t399 - t362 * t416) * r_i_i_C(1) + (-t352 * t402 + t362 * t417) * r_i_i_C(2) - t352 * pkin(8) - t379 * qJ(3) + t384 * qJD(3) + t427 * t380 + t406 * t351 + t405 * t361, t380, 0, t426 * t342 + t410 * t407 + t411 * (-t348 * qJD(5) - t352 * t400), (-t342 * t399 + t351 * t402) * r_i_i_C(1) + (-t342 * t402 - t351 * t399) * r_i_i_C(2) + ((-t348 * t402 - t361 * t399) * r_i_i_C(1) + (t348 * t399 - t361 * t402) * r_i_i_C(2)) * qJD(6); 0, (t373 * t399 - t378 * t416) * r_i_i_C(1) + (t373 * t402 + t378 * t417) * r_i_i_C(2) + t373 * pkin(8) + (t401 * qJD(3) + (qJ(3) * t404 + t427 * t401) * qJD(2)) * t395 - t406 * t374 + t405 * t377, t414, 0, t426 * t346 + t412 * t407 + t411 * (-t368 * qJD(5) + t373 * t400), (-t346 * t399 - t374 * t402) * r_i_i_C(1) + (-t346 * t402 + t374 * t399) * r_i_i_C(2) + ((-t368 * t402 - t377 * t399) * r_i_i_C(1) + (t368 * t399 - t377 * t402) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end