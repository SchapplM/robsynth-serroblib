% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->14), mult. (83->33), div. (0->0), fcn. (72->8), ass. (0->16)
	t89 = sin(pkin(11));
	t92 = cos(pkin(11));
	t95 = sin(qJ(2));
	t96 = cos(qJ(2));
	t98 = t89 * t96 + t92 * t95;
	t88 = t98 * qJD(2);
	t94 = cos(pkin(6));
	t100 = t94 * t95;
	t99 = pkin(2) * qJD(2);
	t97 = t89 * t95 - t92 * t96;
	t87 = t97 * qJD(2);
	t93 = cos(pkin(10));
	t90 = sin(pkin(10));
	t86 = t94 * t88;
	t85 = t94 * t87;
	t1 = [0, (t90 * t86 + t93 * t87) * r_i_i_C(1) + (-t90 * t85 + t93 * t88) * r_i_i_C(2) + (t90 * t100 - t93 * t96) * t99, 0, 0, 0, 0; 0, (-t93 * t86 + t90 * t87) * r_i_i_C(1) + (t93 * t85 + t90 * t88) * r_i_i_C(2) + (-t93 * t100 - t90 * t96) * t99, 0, 0, 0, 0; 0, (-t95 * pkin(2) - t98 * r_i_i_C(1) + t97 * r_i_i_C(2)) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (59->23), mult. (202->41), div. (0->0), fcn. (192->8), ass. (0->23)
	t162 = sin(pkin(11));
	t165 = cos(pkin(11));
	t168 = sin(qJ(2));
	t169 = cos(qJ(2));
	t161 = t168 * t162 - t165 * t169;
	t159 = t161 * qJD(2);
	t175 = pkin(3) - r_i_i_C(2);
	t174 = r_i_i_C(3) + qJ(4);
	t173 = qJD(2) * pkin(2);
	t167 = cos(pkin(6));
	t172 = t167 * t168;
	t171 = t162 * t169 + t165 * t168;
	t158 = t171 * t167;
	t160 = t171 * qJD(2);
	t170 = qJD(2) * t158;
	t166 = cos(pkin(10));
	t164 = sin(pkin(6));
	t163 = sin(pkin(10));
	t157 = t167 * t159;
	t155 = t164 * t160;
	t152 = t159 * t166 + t163 * t170;
	t150 = t163 * t159 - t166 * t170;
	t1 = [0, -(t158 * t163 + t161 * t166) * qJD(4) - t174 * (-t157 * t163 + t160 * t166) + t175 * t152 + (t163 * t172 - t166 * t169) * t173, 0, -t152, 0, 0; 0, -(-t158 * t166 + t161 * t163) * qJD(4) - t174 * (t157 * t166 + t160 * t163) + t175 * t150 + (-t163 * t169 - t166 * t172) * t173, 0, -t150, 0, 0; 0, -t175 * t155 + (t171 * qJD(4) - t174 * t159 - t168 * t173) * t164, 0, t155, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:20
	% EndTime: 2019-10-09 21:26:20
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (142->40), mult. (472->79), div. (0->0), fcn. (472->10), ass. (0->33)
	t243 = sin(qJ(5));
	t245 = cos(qJ(5));
	t253 = t245 * r_i_i_C(1) - t243 * r_i_i_C(2);
	t260 = t253 * qJD(5) + qJD(4);
	t259 = pkin(2) * qJD(2);
	t239 = sin(pkin(6));
	t258 = t239 * t243;
	t257 = t239 * t245;
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t256 = t242 * t244;
	t254 = pkin(3) + pkin(8) + r_i_i_C(3);
	t237 = sin(pkin(11));
	t240 = cos(pkin(11));
	t246 = cos(qJ(2));
	t252 = t246 * t237 + t244 * t240;
	t235 = t244 * t237 - t246 * t240;
	t251 = t243 * r_i_i_C(1) + t245 * r_i_i_C(2) + qJ(4);
	t250 = t252 * t242;
	t234 = t252 * qJD(2);
	t233 = t235 * qJD(2);
	t249 = qJD(2) * t250;
	t248 = t242 * t233;
	t241 = cos(pkin(10));
	t238 = sin(pkin(10));
	t232 = t235 * t242;
	t230 = t235 * t239;
	t228 = t239 * t234;
	t224 = -t238 * t232 + t241 * t252;
	t222 = t241 * t232 + t238 * t252;
	t219 = t241 * t233 + t238 * t249;
	t217 = t238 * t233 - t241 * t249;
	t1 = [0, (t238 * t256 - t241 * t246) * t259 - t260 * (t241 * t235 + t238 * t250) - t251 * (t241 * t234 - t238 * t248) + t254 * t219, 0, -t219, -t253 * t219 + ((-t224 * t243 - t238 * t257) * r_i_i_C(1) + (-t224 * t245 + t238 * t258) * r_i_i_C(2)) * qJD(5), 0; 0, (-t238 * t246 - t241 * t256) * t259 - t260 * (t238 * t235 - t241 * t250) - t251 * (t238 * t234 + t241 * t248) + t254 * t217, 0, -t217, -t253 * t217 + ((-t222 * t243 + t241 * t257) * r_i_i_C(1) + (-t222 * t245 - t241 * t258) * r_i_i_C(2)) * qJD(5), 0; 0, -t254 * t228 + (t260 * t252 + (-t244 * pkin(2) - t251 * t235) * qJD(2)) * t239, 0, t228, t253 * t228 + ((-t230 * t243 - t242 * t245) * r_i_i_C(1) + (-t230 * t245 + t242 * t243) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:21
	% EndTime: 2019-10-09 21:26:22
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (456->79), mult. (1438->146), div. (0->0), fcn. (1546->12), ass. (0->58)
	t382 = cos(pkin(6));
	t377 = sin(pkin(11));
	t380 = cos(pkin(11));
	t385 = sin(qJ(2));
	t388 = cos(qJ(2));
	t400 = t388 * t377 + t385 * t380;
	t407 = t400 * t382;
	t384 = sin(qJ(5));
	t387 = cos(qJ(5));
	t383 = sin(qJ(6));
	t386 = cos(qJ(6));
	t401 = t383 * r_i_i_C(1) + t386 * r_i_i_C(2);
	t393 = qJD(6) * t401;
	t402 = t386 * r_i_i_C(1) - t383 * r_i_i_C(2);
	t399 = pkin(5) + t402;
	t421 = pkin(9) + r_i_i_C(3);
	t427 = qJD(4) - t384 * t393 + (t421 * t384 + t399 * t387) * qJD(5);
	t368 = t385 * t377 - t388 * t380;
	t394 = qJD(6) * t402;
	t366 = t368 * qJD(2);
	t405 = qJD(2) * t388;
	t406 = qJD(2) * t385;
	t367 = -t377 * t405 - t380 * t406;
	t378 = sin(pkin(10));
	t381 = cos(pkin(10));
	t398 = (t377 * t406 - t380 * t405) * t382;
	t345 = t381 * t367 + t378 * t398;
	t353 = -t381 * t368 - t378 * t407;
	t425 = -t378 * t367 + t381 * t398;
	t424 = t378 * t368 - t381 * t407;
	t422 = t399 * t384 - t421 * t387 + qJ(4);
	t420 = pkin(2) * qJD(2);
	t379 = sin(pkin(6));
	t417 = t379 * t384;
	t416 = t379 * t387;
	t413 = t382 * t384;
	t412 = t382 * t385;
	t404 = qJD(5) * t387;
	t403 = t378 * t417;
	t364 = t368 * t379;
	t355 = t364 * t384 + t382 * t387;
	t397 = pkin(3) + pkin(8) + t401;
	t392 = t368 * t382;
	t349 = -t378 * t400 - t381 * t392;
	t396 = t349 * t384 + t381 * t416;
	t395 = -t349 * t387 + t381 * t417;
	t352 = t378 * t392 - t381 * t400;
	t337 = -t352 * t384 + t378 * t416;
	t365 = t400 * t379;
	t391 = qJD(2) * t407;
	t363 = t379 * t366;
	t362 = qJD(2) * t365;
	t344 = t381 * t366 + t378 * t391;
	t341 = t378 * t366 - t381 * t391;
	t334 = qJD(5) * t413 - t362 * t384 - t364 * t404;
	t332 = t395 * qJD(5) - t341 * t384;
	t330 = qJD(5) * t403 + t344 * t384 + t352 * t404;
	t1 = [0, t352 * t394 + t397 * t344 + (t378 * t412 - t381 * t388) * t420 + t422 * t345 + t427 * t353, 0, -t344, -t421 * t330 - (-t352 * t387 - t403) * t393 + t399 * (-t337 * qJD(5) - t344 * t387), (t330 * t383 + t345 * t386) * r_i_i_C(1) + (t330 * t386 - t345 * t383) * r_i_i_C(2) + ((-t337 * t386 - t353 * t383) * r_i_i_C(1) + (t337 * t383 - t353 * t386) * r_i_i_C(2)) * qJD(6); 0, t349 * t394 + t397 * t341 + (-t378 * t388 - t381 * t412) * t420 - t422 * t425 - t427 * t424, 0, -t341, t421 * t332 - t395 * t393 + t399 * (t396 * qJD(5) - t341 * t387), (-t332 * t383 - t386 * t425) * r_i_i_C(1) + (-t332 * t386 + t383 * t425) * r_i_i_C(2) + ((t383 * t424 + t386 * t396) * r_i_i_C(1) + (-t383 * t396 + t386 * t424) * r_i_i_C(2)) * qJD(6); 0, -t379 * pkin(2) * t406 - t397 * t362 - t422 * t363 - t364 * t394 + t427 * t365, 0, t362, -t421 * t334 - (t364 * t387 - t413) * t393 + t399 * (-t355 * qJD(5) + t362 * t387), (t334 * t383 - t363 * t386) * r_i_i_C(1) + (t334 * t386 + t363 * t383) * r_i_i_C(2) + ((-t355 * t386 - t365 * t383) * r_i_i_C(1) + (t355 * t383 - t365 * t386) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end