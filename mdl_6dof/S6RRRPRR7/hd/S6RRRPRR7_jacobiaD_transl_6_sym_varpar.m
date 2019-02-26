% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:21
% EndTime: 2019-02-26 22:19:22
% DurationCPUTime: 0.72s
% Computational Cost: add. (1275->132), mult. (1850->203), div. (0->0), fcn. (1829->14), ass. (0->84)
t502 = pkin(11) + r_i_i_C(3);
t448 = sin(qJ(6));
t452 = cos(qJ(6));
t465 = r_i_i_C(1) * t452 - r_i_i_C(2) * t448;
t463 = pkin(5) + t465;
t451 = sin(qJ(1));
t447 = cos(pkin(6));
t467 = qJD(2) * t447 + qJD(1);
t450 = sin(qJ(2));
t487 = t450 * t451;
t473 = t447 * t487;
t481 = qJD(2) * t450;
t454 = cos(qJ(2));
t455 = cos(qJ(1));
t484 = t454 * t455;
t409 = -qJD(1) * t473 - t451 * t481 + t467 * t484;
t444 = qJD(3) + qJD(5);
t446 = sin(pkin(6));
t488 = t446 * t455;
t511 = t444 * t488 - t409;
t477 = qJD(6) * t452;
t478 = qJD(6) * t448;
t510 = -r_i_i_C(1) * t478 - t477 * r_i_i_C(2);
t485 = t451 * t454;
t486 = t450 * t455;
t422 = t447 * t486 + t485;
t445 = qJ(3) + pkin(12);
t441 = qJ(5) + t445;
t437 = sin(t441);
t438 = cos(t441);
t412 = -t422 * t438 + t437 * t488;
t472 = t447 * t484;
t421 = -t472 + t487;
t509 = -t412 * t448 - t421 * t452;
t508 = t412 * t452 - t421 * t448;
t429 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t445);
t425 = t429 * qJD(3);
t507 = -(t463 * t437 - t502 * t438) * t444 - t425;
t430 = pkin(4) * cos(t445) + cos(qJ(3)) * pkin(3);
t505 = r_i_i_C(1) * t448 + r_i_i_C(2) * t452;
t428 = pkin(2) + t430;
t503 = t502 * t437 + t463 * t438 + t428;
t498 = pkin(8) + t429;
t423 = t447 * t485 + t486;
t408 = t423 * qJD(1) + t422 * qJD(2);
t497 = t408 * t448;
t496 = t408 * t452;
t424 = -t473 + t484;
t493 = t424 * t438;
t492 = t437 * t444;
t491 = t446 * t450;
t490 = t446 * t451;
t489 = t446 * t454;
t483 = qJD(1) * t451;
t482 = qJD(1) * t455;
t480 = qJD(2) * t454;
t479 = qJD(6) * t438;
t475 = t444 * t491;
t471 = t446 * t483;
t470 = t446 * t482;
t469 = t446 * t481;
t468 = t511 * t438;
t407 = t422 * qJD(1) + t423 * qJD(2);
t464 = t444 * t490 - t407;
t462 = -t422 * t444 + t471;
t461 = t444 * t447 + t446 * t480;
t397 = t511 * t437 + t462 * t438;
t405 = -t437 * t475 + t461 * t438;
t459 = t510 * (-t437 * t491 + t438 * t447) + t502 * t405 + t463 * (-t461 * t437 - t438 * t475);
t398 = -t422 * t492 + t437 * t471 - t468;
t458 = t510 * (-t422 * t437 - t438 * t488) + t502 * t398 + t463 * t397;
t395 = t464 * t437 - t438 * t470 + t444 * t493;
t396 = -t424 * t492 + t437 * t470 + t464 * t438;
t457 = t510 * (-t424 * t437 + t438 * t490) + t502 * t396 - t463 * t395;
t456 = t505 * t479 - t507;
t443 = -pkin(10) - qJ(4) - pkin(9);
t426 = t430 * qJD(3);
t418 = t437 * t447 + t438 * t491;
t414 = t437 * t490 + t493;
t406 = -qJD(1) * t472 - t455 * t480 + t467 * t487;
t400 = -t462 * t437 + t468;
t388 = t396 * t452 - t406 * t448 + (-t414 * t448 + t423 * t452) * qJD(6);
t387 = -t396 * t448 - t406 * t452 + (-t414 * t452 - t423 * t448) * qJD(6);
t1 = [(t400 * t452 - t497) * r_i_i_C(1) + (-t400 * t448 - t496) * r_i_i_C(2) + t400 * pkin(5) - t409 * t428 + t422 * t425 + t408 * t443 - t421 * qJD(4) + t426 * t488 + t502 * t397 + (t509 * r_i_i_C(1) - t508 * r_i_i_C(2)) * qJD(6) + (-t455 * pkin(1) - t498 * t490) * qJD(1) (-t407 * t448 + t424 * t477) * r_i_i_C(1) + (-t407 * t452 - t424 * t478) * r_i_i_C(2) + t407 * t443 + t424 * qJD(4) + t503 * t406 + t456 * t423, t407 * t429 - t424 * t426 + (-t425 * t451 + t430 * t482) * t446 + t457, -t406, t457, r_i_i_C(1) * t387 - t388 * r_i_i_C(2); t426 * t490 + t396 * pkin(5) + t388 * r_i_i_C(1) + t387 * r_i_i_C(2) + t423 * qJD(4) + t406 * t443 - t407 * t428 - t424 * t425 + t502 * t395 + (-pkin(1) * t451 + t498 * t488) * qJD(1) (t409 * t448 + t422 * t477) * r_i_i_C(1) + (t409 * t452 - t422 * t478) * r_i_i_C(2) - t409 * t443 + t422 * qJD(4) - t503 * t408 + t456 * t421, -t409 * t429 - t422 * t426 + (t425 * t455 + t430 * t483) * t446 + t458, t408, t458 (-t398 * t448 + t496) * r_i_i_C(1) + (-t398 * t452 - t497) * r_i_i_C(2) + (t508 * r_i_i_C(1) + t509 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t503 + t465 * qJD(6) + qJD(4)) * t450 + (-qJD(2) * t443 + t505 * (qJD(2) - t479) + t507) * t454) * t446, -t425 * t447 + (-t426 * t450 - t429 * t480) * t446 + t459, t469, t459 (-t405 * t448 + t452 * t469) * r_i_i_C(1) + (-t405 * t452 - t448 * t469) * r_i_i_C(2) + ((-t418 * t452 + t448 * t489) * r_i_i_C(1) + (t418 * t448 + t452 * t489) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
