% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:59
% EndTime: 2019-02-26 19:40:59
% DurationCPUTime: 0.49s
% Computational Cost: add. (737->77), mult. (2326->142), div. (0->0), fcn. (2700->14), ass. (0->68)
t470 = sin(pkin(12));
t471 = sin(pkin(11));
t456 = t471 * t470;
t474 = cos(pkin(12));
t475 = cos(pkin(11));
t463 = t475 * t474;
t477 = cos(pkin(6));
t438 = -t477 * t463 + t456;
t476 = cos(pkin(7));
t436 = t438 * t476;
t472 = sin(pkin(7));
t473 = sin(pkin(6));
t461 = t473 * t472;
t447 = t475 * t461;
t482 = t436 + t447;
t458 = t471 * t474;
t462 = t475 * t470;
t439 = t477 * t458 + t462;
t457 = t471 * t473;
t481 = t439 * t476 - t472 * t457;
t464 = t476 * t473;
t480 = t474 * t464 + t472 * t477;
t419 = t477 * t462 + t458;
t430 = sin(qJ(3));
t478 = cos(qJ(3));
t402 = t419 * t430 + t482 * t478;
t479 = t481 * t478;
t460 = t473 * t470;
t410 = t430 * t460 - t480 * t478;
t428 = sin(qJ(6));
t431 = cos(qJ(6));
t466 = t431 * r_i_i_C(1) - t428 * r_i_i_C(2);
t444 = t466 * qJD(6) + qJD(5);
t469 = qJD(3) * t430;
t468 = pkin(4) + pkin(10) + r_i_i_C(3);
t467 = t419 * t478;
t465 = -r_i_i_C(1) * t428 - r_i_i_C(2) * t431;
t403 = -t482 * t430 + t467;
t412 = t438 * t472 - t475 * t464;
t429 = sin(qJ(4));
t432 = cos(qJ(4));
t455 = t403 * t432 + t412 * t429;
t394 = t403 * t429 - t412 * t432;
t420 = -t477 * t456 + t463;
t405 = t420 * t478 - t481 * t430;
t413 = t439 * t472 + t476 * t457;
t454 = t405 * t432 + t413 * t429;
t396 = t405 * t429 - t413 * t432;
t411 = t480 * t430 + t478 * t460;
t418 = -t474 * t461 + t477 * t476;
t453 = t411 * t432 + t418 * t429;
t406 = t411 * t429 - t418 * t432;
t452 = qJ(5) - t465;
t450 = -pkin(5) - pkin(9) - t466;
t449 = qJD(6) * t465;
t440 = -t452 * t429 - t468 * t432 - pkin(3);
t433 = -t444 * t429 + (t468 * t429 - t452 * t432) * qJD(4);
t409 = t411 * qJD(3);
t408 = t410 * qJD(3);
t404 = t420 * t430 + t479;
t401 = t405 * qJD(3);
t400 = t479 * qJD(3) + t420 * t469;
t399 = -t447 * t469 + (-t430 * t436 + t467) * qJD(3);
t398 = t402 * qJD(3);
t392 = t453 * qJD(4) - t408 * t429;
t390 = t454 * qJD(4) - t400 * t429;
t388 = t455 * qJD(4) - t398 * t429;
t1 = [0, 0, t450 * t400 + t440 * t401 + t433 * t404 + t405 * t449, t444 * t454 + t452 * (-t396 * qJD(4) - t400 * t432) - t468 * t390, t390 (t390 * t431 - t401 * t428) * r_i_i_C(1) + (-t390 * t428 - t401 * t431) * r_i_i_C(2) + ((-t396 * t428 - t404 * t431) * r_i_i_C(1) + (-t396 * t431 + t404 * t428) * r_i_i_C(2)) * qJD(6); 0, 0, t450 * t398 + t440 * t399 + t433 * t402 + t403 * t449, t444 * t455 + t452 * (-t394 * qJD(4) - t398 * t432) - t468 * t388, t388 (t388 * t431 - t399 * t428) * r_i_i_C(1) + (-t388 * t428 - t399 * t431) * r_i_i_C(2) + ((-t394 * t428 - t402 * t431) * r_i_i_C(1) + (-t394 * t431 + t402 * t428) * r_i_i_C(2)) * qJD(6); 0, 0, t450 * t408 + t440 * t409 + t433 * t410 + t411 * t449, t444 * t453 + t452 * (-t406 * qJD(4) - t408 * t432) - t468 * t392, t392 (t392 * t431 - t409 * t428) * r_i_i_C(1) + (-t392 * t428 - t409 * t431) * r_i_i_C(2) + ((-t406 * t428 - t410 * t431) * r_i_i_C(1) + (-t406 * t431 + t410 * t428) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
