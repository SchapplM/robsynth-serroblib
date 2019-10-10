% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
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
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
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
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (116->40), mult. (277->78), div. (0->0), fcn. (252->10), ass. (0->33)
	t218 = r_i_i_C(3) + qJ(4) + pkin(8);
	t195 = sin(pkin(11));
	t196 = sin(pkin(6));
	t217 = t195 * t196;
	t197 = cos(pkin(11));
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
	t194 = qJ(3) + pkin(12);
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
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (268->47), mult. (379->81), div. (0->0), fcn. (343->12), ass. (0->43)
	t266 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
	t236 = qJ(3) + pkin(12);
	t233 = qJ(5) + t236;
	t229 = sin(t233);
	t235 = qJD(3) + qJD(5);
	t265 = t229 * t235;
	t230 = cos(t233);
	t264 = t230 * t235;
	t237 = sin(pkin(11));
	t238 = sin(pkin(6));
	t263 = t237 * t238;
	t239 = cos(pkin(11));
	t262 = t238 * t239;
	t240 = cos(pkin(6));
	t242 = sin(qJ(2));
	t261 = t240 * t242;
	t244 = cos(qJ(2));
	t260 = t240 * t244;
	t220 = t237 * t244 + t239 * t261;
	t255 = qJD(2) * t244;
	t252 = t239 * t255;
	t256 = qJD(2) * t242;
	t253 = t237 * t256;
	t215 = -t240 * t252 + t253;
	t250 = t235 * t262 + t215;
	t259 = (-t220 * t264 + t250 * t229) * r_i_i_C(1) + (t220 * t265 + t250 * t230) * r_i_i_C(2);
	t222 = -t237 * t261 + t239 * t244;
	t247 = t237 * t260 + t239 * t242;
	t217 = t247 * qJD(2);
	t249 = -t235 * t263 + t217;
	t258 = (-t222 * t264 + t249 * t229) * r_i_i_C(1) + (t222 * t265 + t249 * t230) * r_i_i_C(2);
	t246 = -t235 * t240 - t238 * t255;
	t254 = t235 * t238 * t242;
	t257 = (t246 * t229 - t230 * t254) * r_i_i_C(1) + (t229 * t254 + t246 * t230) * r_i_i_C(2);
	t226 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t236);
	t251 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t236);
	t248 = -r_i_i_C(1) * t230 + r_i_i_C(2) * t229 - pkin(2) + t251;
	t223 = t226 * qJD(3);
	t245 = t223 + (-r_i_i_C(1) * t229 - r_i_i_C(2) * t230) * t235;
	t224 = t251 * qJD(3);
	t218 = -t240 * t253 + t252;
	t216 = t220 * qJD(2);
	t1 = [0, t222 * qJD(4) - t266 * t217 + t248 * t218 - t245 * t247, -t217 * t226 + t222 * t224 + t223 * t263 + t258, t218, t258, 0; 0, t220 * qJD(4) - t266 * t215 + t245 * (-t237 * t242 + t239 * t260) + t248 * t216, -t215 * t226 + t220 * t224 - t223 * t262 + t259, t216, t259, 0; 0, (qJD(4) * t242 + t245 * t244 + (t248 * t242 + t266 * t244) * qJD(2)) * t238, t240 * t223 + (t224 * t242 + t226 * t255) * t238 + t257, t238 * t256, t257, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:38
	% EndTime: 2019-10-09 22:25:38
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (813->94), mult. (1084->159), div. (0->0), fcn. (1067->14), ass. (0->64)
	t416 = pkin(10) + r_i_i_C(3);
	t373 = sin(qJ(6));
	t376 = cos(qJ(6));
	t389 = t376 * r_i_i_C(1) - t373 * r_i_i_C(2);
	t422 = pkin(5) + t389;
	t401 = qJD(6) * t376;
	t402 = qJD(6) * t373;
	t421 = -r_i_i_C(1) * t402 - t401 * r_i_i_C(2);
	t370 = qJ(3) + pkin(12);
	t359 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t370);
	t356 = t359 * qJD(3);
	t367 = qJ(5) + t370;
	t363 = sin(t367);
	t364 = cos(t367);
	t369 = qJD(3) + qJD(5);
	t420 = -(t363 * t422 - t416 * t364) * t369 + t356;
	t375 = sin(qJ(2));
	t378 = cos(qJ(2));
	t371 = sin(pkin(11));
	t412 = cos(pkin(6));
	t396 = t371 * t412;
	t411 = cos(pkin(11));
	t353 = -t375 * t396 + t411 * t378;
	t418 = t373 * r_i_i_C(1) + t376 * r_i_i_C(2);
	t410 = t363 * t369;
	t409 = t364 * t369;
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
	t395 = t372 * t411;
	t391 = t364 * t395;
	t390 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t370);
	t352 = t411 * t375 + t378 * t396;
	t348 = t352 * qJD(2);
	t388 = t369 * t408 - t348;
	t387 = t412 * t411;
	t385 = t378 * t387;
	t384 = t412 * t369 + t372 * t404;
	t351 = t371 * t378 + t375 * t387;
	t383 = -t416 * t363 - t364 * t422 - pkin(2) + t390;
	t346 = -qJD(2) * t385 + t371 * t405;
	t330 = -t346 * t364 - t351 * t410 - t369 * t391;
	t382 = t421 * (-t351 * t363 - t391) + t416 * t330 + t422 * (-t351 * t409 + (t369 * t395 + t346) * t363);
	t332 = -t353 * t410 + t388 * t364;
	t381 = t421 * (-t353 * t363 + t364 * t408) + t416 * t332 + t422 * (-t353 * t409 - t388 * t363);
	t337 = t384 * t364 - t369 * t399;
	t380 = t421 * (t412 * t364 - t399) + t416 * t337 + t422 * (-t384 * t363 - t369 * t398);
	t379 = t418 * t403 - t420;
	t368 = -pkin(9) - qJ(4) - pkin(8);
	t357 = t390 * qJD(3);
	t350 = t371 * t375 - t385;
	t349 = t353 * qJD(2);
	t347 = t351 * qJD(2);
	t345 = t412 * t363 + t398;
	t341 = t353 * t364 + t363 * t408;
	t339 = t351 * t364 - t363 * t395;
	t1 = [0, (-t348 * t373 + t353 * t401) * r_i_i_C(1) + (-t348 * t376 - t353 * t402) * r_i_i_C(2) + t348 * t368 + t353 * qJD(4) + t383 * t349 + t379 * t352, -t348 * t359 + t353 * t357 + t356 * t408 + t381, t349, t381, (-t332 * t373 + t349 * t376) * r_i_i_C(1) + (-t332 * t376 - t349 * t373) * r_i_i_C(2) + ((-t341 * t376 - t352 * t373) * r_i_i_C(1) + (t341 * t373 - t352 * t376) * r_i_i_C(2)) * qJD(6); 0, (-t346 * t373 + t351 * t401) * r_i_i_C(1) + (-t346 * t376 - t351 * t402) * r_i_i_C(2) + t346 * t368 + t351 * qJD(4) + t383 * t347 + t379 * t350, -t346 * t359 + t351 * t357 - t356 * t395 + t382, t347, t382, (-t330 * t373 + t347 * t376) * r_i_i_C(1) + (-t330 * t376 - t347 * t373) * r_i_i_C(2) + ((-t339 * t376 - t350 * t373) * r_i_i_C(1) + (t339 * t373 - t350 * t376) * r_i_i_C(2)) * qJD(6); 0, ((t383 * qJD(2) + t389 * qJD(6) + qJD(4)) * t375 + (-qJD(2) * t368 + t418 * (qJD(2) - t403) + t420) * t378) * t372, t412 * t356 + (t357 * t375 + t359 * t404) * t372 + t380, t397, t380, (-t337 * t373 + t376 * t397) * r_i_i_C(1) + (-t337 * t376 - t373 * t397) * r_i_i_C(2) + ((-t345 * t376 + t373 * t406) * r_i_i_C(1) + (t345 * t373 + t376 * t406) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end