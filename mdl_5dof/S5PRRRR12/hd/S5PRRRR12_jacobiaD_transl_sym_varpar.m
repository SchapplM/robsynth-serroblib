% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR12
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(5));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(5)) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (119->20), mult. (101->35), div. (0->0), fcn. (60->12), ass. (0->25)
	t119 = qJD(2) + qJD(3);
	t135 = -t119 / 0.2e1;
	t134 = pkin(2) * qJD(2);
	t120 = qJ(2) + qJ(3);
	t116 = pkin(5) - t120;
	t133 = t119 * sin(t116);
	t115 = pkin(5) + t120;
	t132 = t119 * cos(t115);
	t121 = sin(pkin(11));
	t131 = t119 * t121;
	t122 = cos(pkin(11));
	t130 = t119 * t122;
	t124 = sin(qJ(2));
	t129 = cos(pkin(5)) * t124;
	t111 = sin(t115) * t135;
	t109 = t133 / 0.2e1 + t111;
	t112 = cos(t116) * t135;
	t110 = t112 - t132 / 0.2e1;
	t117 = sin(t120);
	t118 = cos(t120);
	t128 = (t122 * t109 - t118 * t131) * r_i_i_C(1) + (t122 * t110 + t117 * t131) * r_i_i_C(2);
	t127 = (-t121 * t109 - t118 * t130) * r_i_i_C(1) + (-t121 * t110 + t117 * t130) * r_i_i_C(2);
	t126 = (t132 / 0.2e1 + t112) * r_i_i_C(1) + (t111 - t133 / 0.2e1) * r_i_i_C(2);
	t125 = cos(qJ(2));
	t1 = [0, (t121 * t129 - t122 * t125) * t134 + t127, t127, 0, 0; 0, (-t121 * t125 - t122 * t129) * t134 + t128, t128, 0, 0; 0, -sin(pkin(5)) * t124 * t134 + t126, t126, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (308->33), mult. (242->51), div. (0->0), fcn. (158->15), ass. (0->38)
	t138 = qJD(2) + qJD(3);
	t135 = qJD(4) + t138;
	t161 = -t135 / 0.2e1;
	t156 = qJ(2) + qJ(3);
	t137 = qJ(4) + t156;
	t131 = pkin(5) - t137;
	t160 = t135 * sin(t131);
	t130 = pkin(5) + t137;
	t159 = t135 * cos(t130);
	t139 = sin(pkin(11));
	t158 = t135 * t139;
	t141 = cos(pkin(11));
	t157 = t135 * t141;
	t126 = sin(t130) * t161;
	t123 = t160 / 0.2e1 + t126;
	t127 = cos(t131) * t161;
	t124 = t127 - t159 / 0.2e1;
	t132 = sin(t137);
	t133 = cos(t137);
	t155 = (t141 * t123 - t133 * t158) * r_i_i_C(1) + (t141 * t124 + t132 * t158) * r_i_i_C(2);
	t154 = (-t139 * t123 - t133 * t157) * r_i_i_C(1) + (-t139 * t124 + t132 * t157) * r_i_i_C(2);
	t153 = (t159 / 0.2e1 + t127) * r_i_i_C(1) + (t126 - t160 / 0.2e1) * r_i_i_C(2);
	t146 = cos(qJ(2));
	t152 = qJD(2) * t146;
	t151 = pkin(3) * t138 * cos(t156);
	t143 = sin(qJ(3));
	t144 = sin(qJ(2));
	t145 = cos(qJ(3));
	t150 = -t143 * t146 - t144 * t145;
	t149 = t150 * qJD(3);
	t148 = pkin(3) * (t150 * qJD(2) + t149);
	t147 = -qJD(2) * (t145 * pkin(3) + pkin(2)) * t144 + (-t143 * t152 + t149) * pkin(3);
	t142 = cos(pkin(5));
	t140 = sin(pkin(5));
	t125 = -pkin(2) * t152 - t151;
	t120 = t142 * t148;
	t119 = t147 * t142;
	t1 = [0, -t139 * t119 + t141 * t125 + t154, -t139 * t120 - t141 * t151 + t154, t154, 0; 0, t141 * t119 + t139 * t125 + t155, t141 * t120 - t139 * t151 + t155, t155, 0; 0, t147 * t140 + t153, t140 * t148 + t153, t153, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:13
	% EndTime: 2024-09-28 18:09:14
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (922->128), mult. (1798->222), div. (0->0), fcn. (1671->16), ass. (0->105)
	t402 = sin(pkin(11));
	t407 = cos(pkin(5));
	t405 = cos(pkin(11));
	t403 = sin(pkin(6));
	t487 = pkin(10) * t403;
	t451 = t405 * t487;
	t386 = -t402 * pkin(4) + t407 * t451;
	t453 = t402 * t487;
	t464 = t405 * t407;
	t388 = pkin(4) * t464 + t453;
	t409 = sin(qJ(4));
	t413 = cos(qJ(4));
	t500 = t386 * t413 - t388 * t409;
	t352 = t500 * qJD(4);
	t410 = sin(qJ(3));
	t411 = sin(qJ(2));
	t415 = cos(qJ(2));
	t401 = qJ(2) + qJ(3) + qJ(4);
	t396 = sin(t401);
	t397 = cos(t401);
	t412 = cos(qJ(5));
	t406 = cos(pkin(6));
	t408 = sin(qJ(5));
	t461 = t407 * t408;
	t450 = t406 * t461;
	t373 = t402 * t412 + t405 * t450;
	t511 = qJD(3) + qJD(2);
	t399 = qJD(4) + t511;
	t462 = t406 * t412;
	t418 = t402 * t462 + t405 * t461;
	t440 = t418 * qJD(5) + t373 * t399;
	t460 = t407 * t412;
	t449 = t406 * t460;
	t374 = -t402 * t408 + t405 * t449;
	t463 = t406 * t408;
	t375 = t402 * t463 - t405 * t460;
	t441 = t375 * qJD(5) - t374 * t399;
	t442 = -t374 * qJD(5) + t375 * t399;
	t443 = -t373 * qJD(5) - t399 * t418;
	t454 = r_i_i_C(3) * t399 * t403;
	t447 = (-t443 * t396 + t441 * t397) * r_i_i_C(2) + (t442 * t396 - t440 * t397) * r_i_i_C(1) + (-t396 * t402 + t397 * t464) * t454;
	t499 = t386 * t409 + t388 * t413;
	t339 = pkin(3) * t464 + t499;
	t344 = -pkin(3) * t402 + t500;
	t414 = cos(qJ(3));
	t433 = -t339 * t410 + t344 * t414;
	t518 = t339 * t414 + t344 * t410;
	t489 = -t518 * t411 + t433 * t415;
	t493 = t499 * qJD(4);
	t491 = t352 * t414 - t410 * t493;
	t502 = t493 * t414;
	t523 = -(t352 * t410 + t502) * t411 + t491 * t415 + t447 + t489 * qJD(3);
	t385 = t405 * pkin(4) + t407 * t453;
	t465 = t402 * t407;
	t387 = pkin(4) * t465 - t451;
	t358 = t385 * t409 + t387 * t413;
	t338 = -pkin(3) * t465 - t358;
	t432 = t385 * t413 - t387 * t409;
	t342 = t405 * pkin(3) + t432;
	t434 = -t338 * t410 - t342 * t414;
	t495 = -t338 * t414 + t342 * t410;
	t490 = t495 * t411 + t434 * t415;
	t350 = t432 * qJD(4);
	t371 = t402 * t450 - t405 * t412;
	t378 = t402 * t461 - t405 * t462;
	t438 = t378 * qJD(5) + t371 * t399;
	t377 = t402 * t460 + t405 * t463;
	t417 = t402 * t449 + t405 * t408;
	t439 = t377 * qJD(5) + t399 * t417;
	t444 = t417 * qJD(5) + t377 * t399;
	t445 = t371 * qJD(5) + t378 * t399;
	t448 = (-t445 * t396 + t439 * t397) * r_i_i_C(2) + (t444 * t396 + t438 * t397) * r_i_i_C(1) + (-t396 * t405 - t397 * t465) * t454;
	t494 = t358 * qJD(4);
	t492 = -t350 * t414 + t410 * t494;
	t501 = t494 * t414;
	t517 = (t350 * t410 + t501) * t411 + t492 * t415 + t448;
	t510 = t358 * t410 - t432 * t414;
	t509 = t358 * t414 + t432 * t410;
	t404 = sin(pkin(5));
	t416 = (r_i_i_C(1) * t408 + r_i_i_C(2) * t412) * qJD(5) * t403;
	t508 = t404 * t416;
	t505 = -t410 * t499 + t500 * t414;
	t504 = t500 * t410 + t414 * t499;
	t394 = -t409 * pkin(4) + t413 * t487;
	t389 = t394 * t414;
	t393 = -pkin(4) * t413 - t409 * t487;
	t392 = pkin(3) - t393;
	t497 = t410 * t392 - t389;
	t390 = t393 * qJD(4);
	t391 = t394 * qJD(4);
	t429 = -t390 * t410 - t391 * t414;
	t484 = (t497 * qJD(3) + t429) * t415;
	t466 = t394 * t410;
	t458 = qJD(2) * t411;
	t457 = qJD(2) * t415;
	t456 = qJD(4) * t410;
	t424 = t396 * t463 - t397 * t412;
	t425 = -t396 * t462 - t397 * t408;
	t426 = -t396 * t412 - t397 * t463;
	t427 = t396 * t408 - t397 * t462;
	t446 = ((t425 * qJD(5) + t426 * t399) * r_i_i_C(1) + (t424 * qJD(5) + t427 * t399) * r_i_i_C(2) + t397 * t454) * t404;
	t437 = t390 * t414 - t410 * t391;
	t428 = -t392 * t414 - t466;
	t435 = t446 + ((t428 * qJD(3) + t437) * t411 - t497 * t457) * t404;
	t1 = [0, t490 * qJD(3) + (-(t405 * pkin(2) - t434) * t415 - (-pkin(2) * t465 - t495) * t411) * qJD(2) + t517, t511 * t490 + t517, (t510 * qJD(3) + t492) * t415 + t509 * t458 + (t509 * qJD(3) + t432 * t456 + t501) * t411 + t510 * t457 + t448, -t402 * t508 + (t445 * r_i_i_C(1) + t444 * r_i_i_C(2)) * t397 + (t439 * r_i_i_C(1) - t438 * r_i_i_C(2)) * t396; 0, (-(t402 * pkin(2) - t433) * t415 - (pkin(2) * t464 + t518) * t411) * qJD(2) + t523, t489 * qJD(2) + t523, (t505 * qJD(3) + t491) * t415 - t504 * t458 + (-t504 * qJD(3) - t500 * t456 - t502) * t411 + t505 * t457 + t447, t405 * t508 + (t443 * r_i_i_C(1) + t442 * r_i_i_C(2)) * t397 + (t441 * r_i_i_C(1) + t440 * r_i_i_C(2)) * t396; 0, (-t484 - (pkin(2) - t428) * t458) * t404 + t435, (t428 * t458 - t484) * t404 + t435, (t437 * t411 - t429 * t415 + t511 * ((t393 * t410 + t389) * t415 + (t393 * t414 - t466) * t411)) * t404 + t446, -t407 * t416 + ((t425 * r_i_i_C(1) + t424 * r_i_i_C(2)) * t399 + (t426 * r_i_i_C(1) + t427 * r_i_i_C(2)) * qJD(5)) * t404;];
	JaD_transl = t1;
end