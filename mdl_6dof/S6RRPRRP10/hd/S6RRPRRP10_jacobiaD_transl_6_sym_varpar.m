% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:07
% EndTime: 2019-02-26 21:51:08
% DurationCPUTime: 0.76s
% Computational Cost: add. (1151->127), mult. (2392->198), div. (0->0), fcn. (2468->12), ass. (0->78)
t471 = cos(pkin(6));
t475 = sin(qJ(1));
t474 = sin(qJ(2));
t517 = t475 * t474;
t507 = t471 * t517;
t512 = qJD(2) * t474;
t477 = cos(qJ(2));
t478 = cos(qJ(1));
t514 = t478 * t477;
t444 = -qJD(1) * t507 - t475 * t512 + (qJD(2) * t471 + qJD(1)) * t514;
t515 = t478 * t474;
t516 = t475 * t477;
t455 = t471 * t515 + t516;
t468 = pkin(11) + qJ(4);
t466 = sin(t468);
t467 = cos(t468);
t470 = sin(pkin(6));
t513 = qJD(1) * t470;
t504 = t475 * t513;
t519 = t470 * t478;
t506 = t467 * t519;
t432 = (-qJD(4) * t455 + t504) * t466 - qJD(4) * t506 + t444 * t467;
t456 = t471 * t516 + t515;
t443 = t456 * qJD(1) + t455 * qJD(2);
t473 = sin(qJ(5));
t476 = cos(qJ(5));
t449 = -t455 * t467 + t466 * t519;
t454 = -t471 * t514 + t517;
t537 = t449 * t476 - t454 * t473;
t542 = t537 * qJD(5) - t432 * t473 + t443 * t476;
t538 = t449 * t473 + t454 * t476;
t541 = t538 * qJD(5) + t432 * t476 + t443 * t473;
t522 = t470 * t474;
t453 = t471 * t466 + t467 * t522;
t518 = t473 * t477;
t536 = -t453 * t476 + t470 * t518;
t465 = cos(pkin(11)) * pkin(3) + pkin(2);
t531 = r_i_i_C(2) + pkin(10);
t535 = t467 * pkin(4) + t531 * t466 + t465;
t534 = -pkin(4) * t466 + t531 * t467;
t510 = qJD(4) * t477;
t533 = (qJD(2) * t467 - qJD(5)) * t474 + t466 * t510;
t529 = r_i_i_C(3) + qJ(6);
t532 = r_i_i_C(1) + pkin(5);
t483 = t529 * t473 + t532 * t476 + pkin(4);
t523 = t467 * t473;
t521 = t470 * t475;
t520 = t470 * t477;
t511 = qJD(4) * t466;
t509 = qJD(5) * t467;
t505 = pkin(3) * sin(pkin(11)) + pkin(8);
t503 = t478 * t513;
t502 = qJD(2) * t520;
t501 = t470 * t512;
t442 = t455 * qJD(1) + t456 * qJD(2);
t496 = t456 * t509 - t442;
t495 = t454 * t509 + t444;
t487 = t507 - t514;
t451 = t466 * t521 - t467 * t487;
t492 = t451 * t476 + t456 * t473;
t491 = -t451 * t473 + t456 * t476;
t490 = (qJD(2) - t509) * t477;
t489 = t466 * t487 + t467 * t521;
t488 = -t466 * t522 + t471 * t467;
t484 = qJD(4) * t534;
t441 = t454 * qJD(1) + t487 * qJD(2);
t482 = -qJD(5) * t487 + t441 * t467 + t456 * t511;
t481 = qJD(5) * t455 - t443 * t467 + t454 * t511;
t480 = qJD(6) * t473 + (-t532 * t473 + t529 * t476) * qJD(5);
t479 = t449 * qJD(4) - t444 * t466 + t467 * t504;
t472 = -pkin(9) - qJ(3);
t446 = t488 * qJD(4) + t467 * t502;
t435 = -t536 * qJD(5) + t446 * t473 - t476 * t501;
t430 = t489 * qJD(4) - t442 * t467 + t466 * t503;
t429 = t451 * qJD(4) - t442 * t466 - t467 * t503;
t420 = t491 * qJD(5) + t430 * t476 - t441 * t473;
t419 = t492 * qJD(5) + t430 * t473 + t441 * t476;
t1 = [t538 * qJD(6) - t432 * pkin(4) - t444 * t465 + t443 * t472 - t454 * qJD(3) + t531 * t479 - t532 * t541 + t529 * t542 + (-t478 * pkin(1) - t505 * t521) * qJD(1) -(t456 * t523 - t476 * t487) * qJD(6) + t442 * t472 - t487 * qJD(3) + t532 * (t496 * t473 + t482 * t476) + t529 * (t482 * t473 - t496 * t476) - t456 * t484 + t535 * t441, -t441, -t483 * t429 + t531 * t430 + t480 * t489, t492 * qJD(6) - t532 * t419 + t529 * t420, t419; -t491 * qJD(6) + t430 * pkin(4) - t442 * t465 + t441 * t472 + t456 * qJD(3) + t531 * t429 + t532 * t420 + t529 * t419 + (-t475 * pkin(1) + t505 * t519) * qJD(1) -(t454 * t523 + t455 * t476) * qJD(6) - t444 * t472 + t455 * qJD(3) + t532 * (t495 * t473 + t481 * t476) + t529 * (t481 * t473 - t495 * t476) - t454 * t484 - t535 * t443, t443, t531 * t432 + t480 * (-t455 * t466 - t506) + t483 * t479, -t537 * qJD(6) + t529 * t541 + t532 * t542, -t542; 0 (t532 * (t473 * t490 - t533 * t476) - t529 * (t533 * t473 + t476 * t490) - (-t467 * t518 + t474 * t476) * qJD(6) + t474 * qJD(3) + t534 * t510 + (-t477 * t472 - t474 * t535) * qJD(2)) * t470, t501, t531 * t446 + t480 * t488 + t483 * (-t453 * qJD(4) - t466 * t502) -t536 * qJD(6) + t529 * (t473 * t501 + t446 * t476 + (-t453 * t473 - t476 * t520) * qJD(5)) - t532 * t435, t435;];
JaD_transl  = t1;
