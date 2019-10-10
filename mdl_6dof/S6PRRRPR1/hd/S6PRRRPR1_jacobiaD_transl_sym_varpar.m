% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:44
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:44
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (305->55), mult. (424->94), div. (0->0), fcn. (381->12), ass. (0->49)
	t236 = qJ(3) + qJ(4);
	t232 = sin(t236);
	t268 = pkin(4) * t232;
	t267 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
	t266 = pkin(3) * qJD(3);
	t231 = pkin(12) + t236;
	t229 = sin(t231);
	t235 = qJD(3) + qJD(4);
	t265 = t229 * t235;
	t230 = cos(t231);
	t264 = t230 * t235;
	t233 = cos(t236);
	t263 = t233 * t235;
	t237 = sin(pkin(11));
	t238 = sin(pkin(6));
	t262 = t237 * t238;
	t239 = cos(pkin(11));
	t261 = t238 * t239;
	t240 = cos(pkin(6));
	t242 = sin(qJ(2));
	t260 = t240 * t242;
	t244 = cos(qJ(2));
	t259 = t240 * t244;
	t220 = t237 * t244 + t239 * t260;
	t254 = qJD(2) * t244;
	t251 = t239 * t254;
	t255 = qJD(2) * t242;
	t252 = t237 * t255;
	t215 = -t240 * t251 + t252;
	t250 = t235 * t261 + t215;
	t258 = (-t220 * t264 + t250 * t229) * r_i_i_C(1) + (t220 * t265 + t250 * t230) * r_i_i_C(2);
	t222 = -t237 * t260 + t239 * t244;
	t247 = t237 * t259 + t239 * t242;
	t217 = t247 * qJD(2);
	t249 = -t235 * t262 + t217;
	t257 = (-t222 * t264 + t249 * t229) * r_i_i_C(1) + (t222 * t265 + t249 * t230) * r_i_i_C(2);
	t246 = -t235 * t240 - t238 * t254;
	t253 = t235 * t238 * t242;
	t256 = (t246 * t229 - t230 * t253) * r_i_i_C(1) + (t229 * t253 + t246 * t230) * r_i_i_C(2);
	t243 = cos(qJ(3));
	t248 = -t243 * pkin(3) - pkin(4) * t233 - r_i_i_C(1) * t230 + r_i_i_C(2) * t229 - pkin(2);
	t241 = sin(qJ(3));
	t223 = -t235 * t268 - t241 * t266;
	t245 = t223 + (-r_i_i_C(1) * t229 - r_i_i_C(2) * t230) * t235;
	t226 = -t241 * pkin(3) - t268;
	t224 = -pkin(4) * t263 - t243 * t266;
	t218 = -t240 * t252 + t251;
	t216 = t220 * qJD(2);
	t1 = [0, t222 * qJD(5) - t267 * t217 + t248 * t218 - t245 * t247, -t217 * t226 + t222 * t224 + t223 * t262 + t257, (-t222 * t263 + t249 * t232) * pkin(4) + t257, t218, 0; 0, t220 * qJD(5) - t267 * t215 + t245 * (-t237 * t242 + t239 * t259) + t248 * t216, -t215 * t226 + t220 * t224 - t223 * t261 + t258, (-t220 * t263 + t250 * t232) * pkin(4) + t258, t216, 0; 0, (qJD(5) * t242 + t245 * t244 + (t248 * t242 + t267 * t244) * qJD(2)) * t238, t240 * t223 + (t224 * t242 + t226 * t254) * t238 + t256, (t246 * t232 - t233 * t253) * pkin(4) + t256, t238 * t255, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:45
	% EndTime: 2019-10-09 22:46:46
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (850->102), mult. (1129->172), div. (0->0), fcn. (1105->14), ass. (0->71)
	t419 = pkin(10) + r_i_i_C(3);
	t373 = sin(qJ(6));
	t376 = cos(qJ(6));
	t390 = t376 * r_i_i_C(1) - t373 * r_i_i_C(2);
	t425 = pkin(5) + t390;
	t401 = qJD(6) * t376;
	t402 = qJD(6) * t373;
	t424 = -r_i_i_C(1) * t402 - t401 * r_i_i_C(2);
	t369 = qJD(3) + qJD(4);
	t374 = sin(qJ(3));
	t414 = pkin(3) * qJD(3);
	t370 = qJ(3) + qJ(4);
	t366 = sin(t370);
	t418 = pkin(4) * t366;
	t355 = -t369 * t418 - t374 * t414;
	t365 = pkin(12) + t370;
	t363 = sin(t365);
	t364 = cos(t365);
	t423 = -(t363 * t425 - t419 * t364) * t369 + t355;
	t375 = sin(qJ(2));
	t378 = cos(qJ(2));
	t371 = sin(pkin(11));
	t413 = cos(pkin(6));
	t396 = t371 * t413;
	t412 = cos(pkin(11));
	t353 = -t375 * t396 + t412 * t378;
	t421 = t373 * r_i_i_C(1) + t376 * r_i_i_C(2);
	t411 = t363 * t369;
	t410 = t364 * t369;
	t367 = cos(t370);
	t409 = t367 * t369;
	t372 = sin(pkin(6));
	t408 = t371 * t372;
	t407 = t372 * t375;
	t406 = t372 * t378;
	t405 = qJD(2) * t375;
	t404 = qJD(2) * t378;
	t403 = qJD(6) * t364;
	t399 = t363 * t407;
	t398 = t364 * t407;
	t397 = t372 * t405;
	t395 = t372 * t412;
	t391 = t364 * t395;
	t352 = t412 * t375 + t378 * t396;
	t348 = t352 * qJD(2);
	t389 = t369 * t408 - t348;
	t388 = t413 * t412;
	t386 = t378 * t388;
	t346 = -qJD(2) * t386 + t371 * t405;
	t385 = t369 * t395 + t346;
	t384 = t413 * t369 + t372 * t404;
	t351 = t371 * t378 + t375 * t388;
	t377 = cos(qJ(3));
	t383 = -t377 * pkin(3) - pkin(4) * t367 - t419 * t363 - t364 * t425 - pkin(2);
	t330 = -t346 * t364 - t351 * t411 - t369 * t391;
	t382 = t424 * (-t351 * t363 - t391) + t419 * t330 + t425 * (-t351 * t410 + t385 * t363);
	t332 = -t353 * t411 + t389 * t364;
	t381 = t424 * (-t353 * t363 + t364 * t408) + t419 * t332 + t425 * (-t353 * t410 - t389 * t363);
	t337 = t384 * t364 - t369 * t399;
	t380 = t424 * (t413 * t364 - t399) + t419 * t337 + t425 * (-t384 * t363 - t369 * t398);
	t379 = t421 * t403 - t423;
	t368 = -qJ(5) - pkin(9) - pkin(8);
	t359 = -t374 * pkin(3) - t418;
	t356 = -pkin(4) * t409 - t377 * t414;
	t350 = t371 * t375 - t386;
	t349 = t353 * qJD(2);
	t347 = t351 * qJD(2);
	t345 = t413 * t363 + t398;
	t341 = t353 * t364 + t363 * t408;
	t339 = t351 * t364 - t363 * t395;
	t1 = [0, (-t348 * t373 + t353 * t401) * r_i_i_C(1) + (-t348 * t376 - t353 * t402) * r_i_i_C(2) + t348 * t368 + t353 * qJD(5) + t383 * t349 + t379 * t352, -t348 * t359 + t353 * t356 + t355 * t408 + t381, (-t353 * t409 - t389 * t366) * pkin(4) + t381, t349, (-t332 * t373 + t349 * t376) * r_i_i_C(1) + (-t332 * t376 - t349 * t373) * r_i_i_C(2) + ((-t341 * t376 - t352 * t373) * r_i_i_C(1) + (t341 * t373 - t352 * t376) * r_i_i_C(2)) * qJD(6); 0, (-t346 * t373 + t351 * t401) * r_i_i_C(1) + (-t346 * t376 - t351 * t402) * r_i_i_C(2) + t346 * t368 + t351 * qJD(5) + t383 * t347 + t379 * t350, -t346 * t359 + t351 * t356 - t355 * t395 + t382, (-t351 * t409 + t385 * t366) * pkin(4) + t382, t347, (-t330 * t373 + t347 * t376) * r_i_i_C(1) + (-t330 * t376 - t347 * t373) * r_i_i_C(2) + ((-t339 * t376 - t350 * t373) * r_i_i_C(1) + (t339 * t373 - t350 * t376) * r_i_i_C(2)) * qJD(6); 0, ((t383 * qJD(2) + t390 * qJD(6) + qJD(5)) * t375 + (-qJD(2) * t368 + t421 * (qJD(2) - t403) + t423) * t378) * t372, t413 * t355 + (t356 * t375 + t359 * t404) * t372 + t380, (-t384 * t366 - t407 * t409) * pkin(4) + t380, t397, (-t337 * t373 + t376 * t397) * r_i_i_C(1) + (-t337 * t376 - t373 * t397) * r_i_i_C(2) + ((-t345 * t376 + t373 * t406) * r_i_i_C(1) + (t345 * t373 + t376 * t406) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end