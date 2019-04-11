% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10V2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:35
% EndTime: 2019-04-11 14:56:36
% DurationCPUTime: 1.22s
% Computational Cost: add. (1406->167), mult. (2206->278), div. (0->0), fcn. (2177->12), ass. (0->122)
t524 = qJ(2) + qJ(3);
t521 = sin(t524);
t526 = sin(qJ(5));
t531 = cos(qJ(5));
t523 = qJD(2) + qJD(3);
t532 = cos(qJ(4));
t566 = qJD(5) * t532 - t523;
t527 = sin(qJ(4));
t592 = qJD(4) * t527;
t621 = -t526 * t592 + t566 * t531;
t634 = t621 * t521;
t622 = t566 * t526 + t531 * t592;
t633 = t622 * t521;
t522 = cos(t524);
t534 = cos(qJ(1));
t596 = t534 * t527;
t529 = sin(qJ(1));
t599 = t529 * t532;
t501 = t522 * t599 - t596;
t603 = t523 * t529;
t578 = t522 * t603;
t593 = qJD(1) * t534;
t550 = -t521 * t593 - t578;
t544 = -qJD(5) * t501 - t550;
t595 = t534 * t532;
t576 = t522 * t595;
t600 = t529 * t527;
t503 = t576 + t600;
t590 = qJD(4) * t534;
t568 = t532 * t590;
t569 = t529 * t592;
t579 = t521 * t603;
t489 = t503 * qJD(1) - t522 * t569 - t532 * t579 - t568;
t608 = t521 * t529;
t561 = qJD(5) * t608 + t489;
t472 = t544 * t526 + t561 * t531;
t591 = qJD(4) * t532;
t548 = t527 * t593 + t529 * t591;
t594 = qJD(1) * t529;
t549 = t527 * t590 + t532 * t594;
t488 = t548 * t522 - t527 * t579 - t549;
t525 = sin(qJ(6));
t530 = cos(qJ(6));
t632 = t472 * t525 - t488 * t530;
t631 = -t472 * t530 - t488 * t525;
t620 = pkin(6) + r_i_i_C(3);
t567 = t523 * t532 - qJD(5);
t558 = t567 * t531;
t587 = qJD(6) * t527;
t570 = t521 * t587;
t630 = t522 * t558 + t570 - t633;
t562 = t530 * r_i_i_C(1) - t525 * r_i_i_C(2);
t583 = t620 * t526;
t629 = t562 * t531 + t583;
t491 = t501 * t531 + t526 * t608;
t500 = t522 * t600 + t595;
t628 = t491 * t525 - t500 * t530;
t627 = t491 * t530 + t500 * t525;
t625 = t525 * r_i_i_C(1) + t530 * r_i_i_C(2);
t598 = t531 * t532;
t606 = t522 * t526;
t498 = t521 * t598 - t606;
t624 = qJD(6) * t498;
t528 = sin(qJ(2));
t614 = pkin(2) * qJD(2);
t584 = t528 * t614;
t617 = pkin(5) * t522;
t618 = pkin(3) * t521;
t623 = (-t617 + t618) * t523 + t584;
t619 = pkin(2) * t528;
t609 = t521 * t526;
t607 = t522 * t523;
t605 = t522 * t531;
t604 = t523 * t527;
t602 = t523 * t534;
t601 = t527 * t531;
t597 = t531 * t534;
t589 = qJD(5) * t534;
t588 = qJD(6) * t525;
t586 = qJD(6) * t530;
t585 = qJD(6) * t531;
t582 = t521 * t602;
t581 = t522 * t602;
t580 = t523 * t597;
t577 = t522 * t596;
t575 = t521 * t594;
t574 = t531 * t594;
t571 = t521 * t589;
t564 = -pkin(3) * t522 - pkin(5) * t521;
t560 = -t521 * t558 + (-t622 + t587) * t522;
t485 = t567 * t605 - t633;
t559 = -t485 - t570;
t557 = t567 * t526;
t556 = -t498 * t594 + t630 * t534;
t496 = t498 * t534;
t555 = qJD(1) * t496 + t630 * t529;
t494 = t503 * t531 + t534 * t609;
t554 = t532 * t609 + t605;
t537 = -t523 * t577 + qJD(6) * t496 + (t527 * t594 - t568) * t521;
t553 = (t556 * t525 + t537 * t530) * r_i_i_C(2) + (t537 * t525 - t556 * t530) * r_i_i_C(1) + pkin(3) * t575 + t620 * (-t523 * t526 * t576 - t571 * t598 + (t526 * t589 + t574) * t522 + (t549 * t526 + t580) * t521);
t538 = -t548 * t521 - t527 * t578 + t529 * t624;
t552 = (t555 * t525 + t538 * t530) * r_i_i_C(2) + (t538 * t525 - t555 * t530) * r_i_i_C(1) + t593 * t617 - t620 * (t554 * t593 + (t522 * t557 + t634) * t529);
t533 = cos(qJ(2));
t551 = qJD(1) * (-t533 * pkin(2) - pkin(1) + t564);
t542 = -qJD(6) * (t522 * t598 + t609) - t521 * t604 + t522 * t591;
t545 = -t523 * t618 + (-t560 * t525 + t542 * t530) * r_i_i_C(2) + (t542 * t525 + t560 * t530) * r_i_i_C(1) + pkin(5) * t607 - t620 * (t521 * t557 - t522 * t621);
t543 = t521 * t591 + t522 * t604 - t624;
t541 = t564 * t523 - t533 * t614;
t539 = (t562 * t526 - t620 * t531) * qJD(5);
t471 = -t561 * t526 + t544 * t531;
t535 = t625 * t585 + t539;
t502 = t577 - t599;
t493 = -t503 * t526 + t521 * t597;
t490 = -t501 * t526 + t531 * t608;
t487 = (-qJD(4) * t522 + qJD(1)) * t596 + (-t582 + (-qJD(1) * t522 + qJD(4)) * t529) * t532;
t486 = t500 * qJD(1) - t522 * t568 + t527 * t582 - t569;
t484 = -t567 * t606 - t634;
t470 = (t487 + t571) * t531 + (-qJD(5) * t503 - t575 + t581) * t526;
t469 = t494 * qJD(5) + t487 * t526 + t521 * t574 - t522 * t580;
t462 = t470 * t530 - t486 * t525 + (-t494 * t525 + t502 * t530) * qJD(6);
t461 = -t470 * t525 - t486 * t530 + (-t494 * t530 - t502 * t525) * qJD(6);
t1 = [t631 * r_i_i_C(1) + t632 * r_i_i_C(2) + t620 * t471 + (t628 * r_i_i_C(1) + t627 * r_i_i_C(2)) * qJD(6) + t623 * t529 + t534 * t551 (-t617 + t619) * t594 + t541 * t534 + t553, -pkin(3) * t581 + (-t522 * t594 - t582) * pkin(5) + t553 (t487 * t525 + t503 * t586) * r_i_i_C(1) + (t487 * t530 - t503 * t588) * r_i_i_C(2) + t629 * t486 + t535 * t502, t620 * t470 + (t469 * t525 - t493 * t586) * r_i_i_C(2) + (-t469 * t530 - t493 * t588) * r_i_i_C(1), t461 * r_i_i_C(1) - t462 * r_i_i_C(2); t462 * r_i_i_C(1) + t461 * r_i_i_C(2) + t620 * t469 + t529 * t551 - t623 * t534 (-t618 - t619) * t593 + t541 * t529 + t552, t550 * pkin(3) - pkin(5) * t579 + t552 (t489 * t525 + t501 * t586) * r_i_i_C(1) + (t489 * t530 - t501 * t588) * r_i_i_C(2) - t629 * t488 + t535 * t500, t620 * t472 + (-t471 * t525 - t490 * t586) * r_i_i_C(2) + (t471 * t530 - t490 * t588) * r_i_i_C(1), -t632 * r_i_i_C(1) + t631 * r_i_i_C(2) + (-t627 * r_i_i_C(1) + t628 * r_i_i_C(2)) * qJD(6); 0, t545 - t584, t545 ((t525 * t532 - t530 * t601) * r_i_i_C(1) + (t525 * t601 + t530 * t532) * r_i_i_C(2) - t527 * t583) * t607 + ((-qJD(4) * t629 + t562 * qJD(6)) * t532 + (t539 + t625 * (-qJD(4) + t585)) * t527) * t521, t620 * t485 + (-t484 * t525 + t554 * t586) * r_i_i_C(2) + (t484 * t530 + t554 * t588) * r_i_i_C(1) (t543 * r_i_i_C(1) + t559 * r_i_i_C(2)) * t530 + (t559 * r_i_i_C(1) - t543 * r_i_i_C(2)) * t525;];
JaD_transl  = t1;
