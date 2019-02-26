% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:36
% EndTime: 2019-02-26 19:41:36
% DurationCPUTime: 0.58s
% Computational Cost: add. (823->95), mult. (2546->168), div. (0->0), fcn. (2964->14), ass. (0->76)
t481 = sin(pkin(12));
t482 = sin(pkin(11));
t468 = t482 * t481;
t485 = cos(pkin(12));
t486 = cos(pkin(11));
t474 = t486 * t485;
t488 = cos(pkin(6));
t443 = -t474 * t488 + t468;
t487 = cos(pkin(7));
t441 = t443 * t487;
t483 = sin(pkin(7));
t484 = sin(pkin(6));
t472 = t484 * t483;
t452 = t486 * t472;
t497 = t441 + t452;
t470 = t482 * t485;
t473 = t486 * t481;
t444 = t470 * t488 + t473;
t469 = t482 * t484;
t496 = t444 * t487 - t469 * t483;
t475 = t487 * t484;
t495 = t475 * t485 + t483 * t488;
t492 = pkin(5) + r_i_i_C(1);
t433 = sin(qJ(5));
t436 = cos(qJ(5));
t449 = qJD(5) * (t436 * r_i_i_C(2) + t433 * t492);
t422 = t473 * t488 + t470;
t435 = sin(qJ(3));
t491 = cos(qJ(3));
t405 = t422 * t435 + t491 * t497;
t494 = t496 * t491;
t471 = t484 * t481;
t413 = t435 * t471 - t491 * t495;
t489 = r_i_i_C(3) + qJ(6) + pkin(10);
t480 = qJD(3) * t435;
t479 = qJD(5) * t433;
t478 = qJD(5) * t436;
t477 = t422 * t491;
t401 = t405 * qJD(3);
t437 = cos(qJ(4));
t406 = -t435 * t497 + t477;
t415 = t443 * t483 - t475 * t486;
t434 = sin(qJ(4));
t462 = -t406 * t434 + t415 * t437;
t392 = qJD(4) * t462 - t401 * t437;
t402 = -t452 * t480 + (-t435 * t441 + t477) * qJD(3);
t467 = -t392 * t433 + t402 * t436;
t423 = -t468 * t488 + t474;
t403 = qJD(3) * t494 + t423 * t480;
t408 = t423 * t491 - t435 * t496;
t416 = t444 * t483 + t469 * t487;
t461 = -t408 * t434 + t416 * t437;
t394 = qJD(4) * t461 - t403 * t437;
t404 = t408 * qJD(3);
t466 = -t394 * t433 + t404 * t436;
t411 = t413 * qJD(3);
t414 = t435 * t495 + t471 * t491;
t421 = -t472 * t485 + t487 * t488;
t459 = -t414 * t434 + t421 * t437;
t396 = qJD(4) * t459 - t411 * t437;
t412 = t414 * qJD(3);
t465 = -t396 * t433 + t412 * t436;
t398 = t406 * t437 + t415 * t434;
t464 = -t398 * t436 - t405 * t433;
t400 = t408 * t437 + t416 * t434;
t407 = t423 * t435 + t494;
t463 = -t400 * t436 - t407 * t433;
t410 = t414 * t437 + t421 * t434;
t460 = -t410 * t436 - t413 * t433;
t458 = -t433 * r_i_i_C(2) + t436 * t492 + pkin(4);
t446 = -t434 * t489 - t437 * t458 - pkin(3);
t438 = -t434 * qJD(6) + t437 * t449 + (t434 * t458 - t437 * t489) * qJD(4);
t395 = qJD(4) * t410 - t411 * t434;
t393 = qJD(4) * t400 - t403 * t434;
t391 = qJD(4) * t398 - t401 * t434;
t1 = [0, 0 (-t403 * t436 - t408 * t479) * r_i_i_C(2) - t403 * pkin(9) + t446 * t404 + t438 * t407 + t492 * (-t403 * t433 + t408 * t478) t400 * qJD(6) - t393 * t458 + t394 * t489 - t449 * t461, t466 * r_i_i_C(1) + (-t394 * t436 - t404 * t433) * r_i_i_C(2) + (t463 * r_i_i_C(1) + (t400 * t433 - t407 * t436) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t463 + t466) * pkin(5), t393; 0, 0 (-t401 * t436 - t406 * t479) * r_i_i_C(2) - t401 * pkin(9) + t446 * t402 + t438 * t405 + t492 * (-t401 * t433 + t406 * t478) t398 * qJD(6) - t391 * t458 + t392 * t489 - t449 * t462, t467 * r_i_i_C(1) + (-t392 * t436 - t402 * t433) * r_i_i_C(2) + (t464 * r_i_i_C(1) + (t398 * t433 - t405 * t436) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t464 + t467) * pkin(5), t391; 0, 0 (-t411 * t436 - t414 * t479) * r_i_i_C(2) - t411 * pkin(9) + t446 * t412 + t438 * t413 + t492 * (-t411 * t433 + t414 * t478) t410 * qJD(6) - t395 * t458 + t396 * t489 - t449 * t459, t465 * r_i_i_C(1) + (-t396 * t436 - t412 * t433) * r_i_i_C(2) + (t460 * r_i_i_C(1) + (t410 * t433 - t413 * t436) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t460 + t465) * pkin(5), t395;];
JaD_transl  = t1;
