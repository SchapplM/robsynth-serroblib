% Calculate time derivative of joint inertia matrix for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:31
% EndTime: 2019-03-09 16:22:36
% DurationCPUTime: 35.58s
% Computational Cost: add. (80053->1480), mult. (197474->1946), div. (0->0), fcn. (219105->12), ass. (0->598)
t555 = cos(pkin(6));
t557 = sin(qJ(2));
t553 = sin(pkin(6));
t732 = cos(qJ(3));
t658 = t553 * t732;
t731 = sin(qJ(3));
t525 = t555 * t731 + t557 * t658;
t733 = cos(qJ(2));
t645 = qJD(2) * t733;
t617 = t553 * t645;
t493 = qJD(3) * t525 + t617 * t731;
t552 = sin(pkin(11));
t554 = cos(pkin(11));
t684 = qJD(2) * t557;
t648 = t553 * t684;
t436 = t493 * t554 - t552 * t648;
t724 = t493 * t552;
t437 = t554 * t648 + t724;
t643 = qJD(3) * t731;
t616 = t553 * t643;
t644 = qJD(3) * t732;
t494 = t555 * t644 - t557 * t616 + t617 * t732;
t309 = Icges(6,5) * t437 + Icges(6,6) * t436 + Icges(6,3) * t494;
t310 = Icges(6,4) * t437 + Icges(6,2) * t436 + Icges(6,6) * t494;
t311 = Icges(6,1) * t437 + Icges(6,4) * t436 + Icges(6,5) * t494;
t657 = t553 * t731;
t524 = -t555 * t732 + t557 * t657;
t659 = t553 * t733;
t491 = t524 * t554 + t552 * t659;
t721 = t524 * t552;
t492 = -t554 * t659 + t721;
t354 = Icges(6,5) * t492 + Icges(6,6) * t491 + Icges(6,3) * t525;
t355 = Icges(6,4) * t492 + Icges(6,2) * t491 + Icges(6,6) * t525;
t356 = Icges(6,1) * t492 + Icges(6,4) * t491 + Icges(6,5) * t525;
t105 = t525 * t309 + t491 * t310 + t492 * t311 + t494 * t354 + t436 * t355 + t437 * t356;
t392 = Icges(5,5) * t648 - Icges(5,6) * t494 + Icges(5,3) * t493;
t393 = Icges(5,4) * t648 - Icges(5,2) * t494 + Icges(5,6) * t493;
t394 = Icges(5,1) * t648 - Icges(5,4) * t494 + Icges(5,5) * t493;
t438 = -Icges(5,5) * t659 - Icges(5,6) * t525 + Icges(5,3) * t524;
t439 = -Icges(5,4) * t659 - Icges(5,2) * t525 + Icges(5,6) * t524;
t440 = -Icges(5,1) * t659 - Icges(5,4) * t525 + Icges(5,5) * t524;
t154 = t524 * t392 - t525 * t393 - t394 * t659 + t493 * t438 - t494 * t439 + t440 * t648;
t395 = Icges(4,5) * t494 - Icges(4,6) * t493 + Icges(4,3) * t648;
t396 = Icges(4,4) * t494 - Icges(4,2) * t493 + Icges(4,6) * t648;
t397 = Icges(4,1) * t494 - Icges(4,4) * t493 + Icges(4,5) * t648;
t441 = Icges(4,5) * t525 - Icges(4,6) * t524 - Icges(4,3) * t659;
t442 = Icges(4,4) * t525 - Icges(4,2) * t524 - Icges(4,6) * t659;
t443 = Icges(4,1) * t525 - Icges(4,4) * t524 - Icges(4,5) * t659;
t155 = -t395 * t659 - t524 * t396 + t525 * t397 + t441 * t648 - t493 * t442 + t494 * t443;
t753 = t105 + t154 + t155;
t558 = sin(qJ(1));
t656 = t558 * t733;
t559 = cos(qJ(1));
t717 = t559 * t557;
t577 = -t555 * t717 - t656;
t620 = t559 * t658;
t495 = -t577 * t731 + t620;
t655 = t559 * t733;
t619 = t555 * t655;
t718 = t558 * t557;
t576 = t619 - t718;
t427 = t495 * t554 + t552 * t576;
t723 = t495 * t552;
t428 = -t554 * t576 + t723;
t654 = t559 * t731;
t496 = -t553 * t654 - t577 * t732;
t288 = Icges(6,5) * t428 + Icges(6,6) * t427 + Icges(6,3) * t496;
t290 = Icges(6,4) * t428 + Icges(6,2) * t427 + Icges(6,6) * t496;
t292 = Icges(6,1) * t428 + Icges(6,4) * t427 + Icges(6,5) * t496;
t132 = t288 * t525 + t290 * t491 + t292 * t492;
t364 = -Icges(5,5) * t576 - Icges(5,6) * t496 + Icges(5,3) * t495;
t368 = -Icges(5,4) * t576 - Icges(5,2) * t496 + Icges(5,6) * t495;
t372 = -Icges(5,1) * t576 - Icges(5,4) * t496 + Icges(5,5) * t495;
t198 = t524 * t364 - t525 * t368 - t372 * t659;
t366 = Icges(4,5) * t496 - Icges(4,6) * t495 - Icges(4,3) * t576;
t370 = Icges(4,4) * t496 - Icges(4,2) * t495 - Icges(4,6) * t576;
t374 = Icges(4,1) * t496 - Icges(4,4) * t495 - Icges(4,5) * t576;
t200 = -t366 * t659 - t524 * t370 + t525 * t374;
t752 = t132 + t198 + t200;
t529 = -t555 * t718 + t655;
t497 = t529 * t731 - t558 * t658;
t578 = -t555 * t656 - t717;
t429 = t497 * t554 + t552 * t578;
t722 = t497 * t552;
t430 = -t554 * t578 + t722;
t621 = t558 * t657;
t498 = t529 * t732 + t621;
t289 = Icges(6,5) * t430 + Icges(6,6) * t429 + Icges(6,3) * t498;
t291 = Icges(6,4) * t430 + Icges(6,2) * t429 + Icges(6,6) * t498;
t293 = Icges(6,1) * t430 + Icges(6,4) * t429 + Icges(6,5) * t498;
t133 = t289 * t525 + t291 * t491 + t293 * t492;
t365 = -Icges(5,5) * t578 - Icges(5,6) * t498 + Icges(5,3) * t497;
t369 = -Icges(5,4) * t578 - Icges(5,2) * t498 + Icges(5,6) * t497;
t373 = -Icges(5,1) * t578 - Icges(5,4) * t498 + Icges(5,5) * t497;
t199 = t524 * t365 - t525 * t369 - t373 * t659;
t367 = Icges(4,5) * t498 - Icges(4,6) * t497 - Icges(4,3) * t578;
t371 = Icges(4,4) * t498 - Icges(4,2) * t497 - Icges(4,6) * t578;
t375 = Icges(4,1) * t498 - Icges(4,4) * t497 - Icges(4,5) * t578;
t201 = -t367 * t659 - t524 * t371 + t525 * t375;
t751 = t133 + t199 + t201;
t469 = qJD(1) * t577 + qJD(2) * t578;
t618 = qJD(1) * t658;
t357 = qJD(3) * t498 + t469 * t731 - t559 * t618;
t551 = pkin(11) + qJ(6);
t548 = sin(t551);
t549 = cos(t551);
t416 = t497 * t548 - t549 * t578;
t468 = -qJD(1) * t619 - t559 * t645 + (qJD(2) * t555 + qJD(1)) * t718;
t248 = -qJD(6) * t416 + t357 * t549 + t468 * t548;
t415 = t497 * t549 + t548 * t578;
t249 = qJD(6) * t415 + t357 * t548 - t468 * t549;
t358 = t469 * t732 - t529 * t643 + (qJD(1) * t654 + t558 * t644) * t553;
t148 = t249 * rSges(7,1) + t248 * rSges(7,2) + t358 * rSges(7,3);
t547 = pkin(5) * t554 + pkin(4);
t556 = -pkin(10) - qJ(5);
t726 = t357 * t552;
t750 = pkin(5) * t726 - t358 * t556 - t468 * t547 + t148;
t282 = t416 * rSges(7,1) + t415 * rSges(7,2) + t498 * rSges(7,3);
t749 = pkin(5) * t722 - t498 * t556 - t547 * t578 + t282;
t748 = t753 * t555;
t413 = t495 * t549 + t548 * t576;
t414 = t495 * t548 - t549 * t576;
t273 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t496;
t275 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t496;
t277 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t496;
t121 = t273 * t496 + t275 * t413 + t277 * t414;
t274 = Icges(7,5) * t416 + Icges(7,6) * t415 + Icges(7,3) * t498;
t276 = Icges(7,4) * t416 + Icges(7,2) * t415 + Icges(7,6) * t498;
t278 = Icges(7,1) * t416 + Icges(7,4) * t415 + Icges(7,5) * t498;
t122 = t274 * t496 + t276 * t413 + t278 * t414;
t475 = t524 * t549 + t548 * t659;
t579 = -t524 * t548 + t549 * t659;
t322 = -Icges(7,5) * t579 + Icges(7,6) * t475 + Icges(7,3) * t525;
t323 = -Icges(7,4) * t579 + Icges(7,2) * t475 + Icges(7,6) * t525;
t324 = -Icges(7,1) * t579 + Icges(7,4) * t475 + Icges(7,5) * t525;
t156 = t322 * t496 + t323 * t413 + t324 * t414;
t471 = qJD(1) * t529 + qJD(2) * t576;
t359 = t471 * t731 - t558 * t618 - t559 * t616 - t577 * t644;
t470 = -qJD(1) * t578 - qJD(2) * t577;
t250 = -qJD(6) * t414 + t359 * t549 - t470 * t548;
t251 = qJD(6) * t413 + t359 * t548 + t470 * t549;
t360 = qJD(1) * t621 - qJD(3) * t620 + t471 * t732 + t577 * t643;
t143 = Icges(7,5) * t251 + Icges(7,6) * t250 + Icges(7,3) * t360;
t145 = Icges(7,4) * t251 + Icges(7,2) * t250 + Icges(7,6) * t360;
t147 = Icges(7,1) * t251 + Icges(7,4) * t250 + Icges(7,5) * t360;
t34 = t143 * t496 + t145 * t413 + t147 * t414 + t250 * t275 + t251 * t277 + t273 * t360;
t142 = Icges(7,5) * t249 + Icges(7,6) * t248 + Icges(7,3) * t358;
t144 = Icges(7,4) * t249 + Icges(7,2) * t248 + Icges(7,6) * t358;
t146 = Icges(7,1) * t249 + Icges(7,4) * t248 + Icges(7,5) * t358;
t35 = t142 * t496 + t144 * t413 + t146 * t414 + t250 * t276 + t251 * t278 + t274 * t360;
t341 = qJD(6) * t579 + t493 * t549 - t548 * t648;
t342 = qJD(6) * t475 + t493 * t548 + t549 * t648;
t236 = Icges(7,5) * t342 + Icges(7,6) * t341 + Icges(7,3) * t494;
t237 = Icges(7,4) * t342 + Icges(7,2) * t341 + Icges(7,6) * t494;
t238 = Icges(7,1) * t342 + Icges(7,4) * t341 + Icges(7,5) * t494;
t62 = t236 * t496 + t237 * t413 + t238 * t414 + t250 * t323 + t251 * t324 + t322 * t360;
t2 = t121 * t360 + t122 * t358 + t156 * t494 + t34 * t496 + t35 * t498 + t525 * t62;
t747 = -t2 / 0.2e1;
t746 = pkin(4) - t547;
t631 = -t468 * pkin(4) + qJ(5) * t358;
t715 = -t631 + t750;
t590 = -t251 * rSges(7,1) - t250 * rSges(7,2);
t149 = t360 * rSges(7,3) - t590;
t344 = t360 * qJ(5);
t725 = t359 * t552;
t677 = pkin(5) * t725;
t714 = -t360 * t556 - t470 * t746 + t149 - t344 + t677;
t589 = -rSges(7,1) * t414 - rSges(7,2) * t413;
t281 = rSges(7,3) * t496 - t589;
t483 = t496 * qJ(5);
t678 = pkin(5) * t723;
t708 = -t496 * t556 + t576 * t746 + t281 - t483 + t678;
t435 = -pkin(4) * t578 + qJ(5) * t498;
t707 = -t435 + t749;
t745 = t470 * pkin(4) + t344;
t550 = t559 * pkin(1);
t720 = t553 * t558;
t688 = pkin(8) * t720 + t550;
t130 = t273 * t525 + t275 * t475 - t277 * t579;
t131 = t274 * t525 + t276 * t475 - t278 * t579;
t38 = t143 * t525 + t145 * t475 - t147 * t579 + t273 * t494 + t275 * t341 + t277 * t342;
t39 = t142 * t525 + t144 * t475 - t146 * t579 + t274 * t494 + t276 * t341 + t278 * t342;
t72 = t525 * t236 + t475 * t237 - t238 * t579 + t494 * t322 + t341 * t323 + t342 * t324;
t71 = t72 * t555;
t11 = t71 + (-t38 * t559 + t39 * t558 + (t130 * t558 + t131 * t559) * qJD(1)) * t553;
t305 = t359 * t554 - t470 * t552;
t306 = t470 * t554 + t725;
t171 = Icges(6,5) * t306 + Icges(6,6) * t305 + Icges(6,3) * t360;
t173 = Icges(6,4) * t306 + Icges(6,2) * t305 + Icges(6,6) * t360;
t175 = Icges(6,1) * t306 + Icges(6,4) * t305 + Icges(6,5) * t360;
t48 = t171 * t525 + t173 * t491 + t175 * t492 + t288 * t494 + t290 * t436 + t292 * t437;
t303 = t357 * t554 + t468 * t552;
t304 = -t468 * t554 + t726;
t170 = Icges(6,5) * t304 + Icges(6,6) * t303 + Icges(6,3) * t358;
t172 = Icges(6,4) * t304 + Icges(6,2) * t303 + Icges(6,6) * t358;
t174 = Icges(6,1) * t304 + Icges(6,4) * t303 + Icges(6,5) * t358;
t49 = t170 * t525 + t172 * t491 + t174 * t492 + t289 * t494 + t291 * t436 + t293 * t437;
t218 = Icges(5,5) * t470 - Icges(5,6) * t360 + Icges(5,3) * t359;
t222 = Icges(5,4) * t470 - Icges(5,2) * t360 + Icges(5,6) * t359;
t226 = Icges(5,1) * t470 - Icges(5,4) * t360 + Icges(5,5) * t359;
t87 = t524 * t218 - t525 * t222 + t493 * t364 - t494 * t368 + (-t226 * t733 + t372 * t684) * t553;
t217 = -Icges(5,5) * t468 - Icges(5,6) * t358 + Icges(5,3) * t357;
t221 = -Icges(5,4) * t468 - Icges(5,2) * t358 + Icges(5,6) * t357;
t225 = -Icges(5,1) * t468 - Icges(5,4) * t358 + Icges(5,5) * t357;
t88 = t524 * t217 - t525 * t221 + t493 * t365 - t494 * t369 + (-t225 * t733 + t373 * t684) * t553;
t220 = Icges(4,5) * t360 - Icges(4,6) * t359 + Icges(4,3) * t470;
t224 = Icges(4,4) * t360 - Icges(4,2) * t359 + Icges(4,6) * t470;
t228 = Icges(4,1) * t360 - Icges(4,4) * t359 + Icges(4,5) * t470;
t89 = -t524 * t224 + t525 * t228 - t493 * t370 + t494 * t374 + (-t220 * t733 + t366 * t684) * t553;
t219 = Icges(4,5) * t358 - Icges(4,6) * t357 - Icges(4,3) * t468;
t223 = Icges(4,4) * t358 - Icges(4,2) * t357 - Icges(4,6) * t468;
t227 = Icges(4,1) * t358 - Icges(4,4) * t357 - Icges(4,5) * t468;
t90 = -t524 * t223 + t525 * t227 - t493 * t371 + t494 * t375 + (-t219 * t733 + t367 * t684) * t553;
t744 = t11 + t748 + ((-t48 - t87 - t89) * t559 + (t49 + t88 + t90) * t558 + (t558 * t752 + t751 * t559) * qJD(1)) * t553;
t743 = -t72 - t753;
t742 = t358 / 0.2e1;
t741 = t360 / 0.2e1;
t740 = t494 / 0.2e1;
t739 = t496 / 0.2e1;
t738 = t498 / 0.2e1;
t737 = t525 / 0.2e1;
t736 = t558 / 0.2e1;
t735 = rSges(5,2) - pkin(3);
t734 = -rSges(6,3) - pkin(3);
t165 = t322 * t525 + t323 * t475 - t324 * t579;
t728 = t165 * t494 + t72 * t525;
t727 = Icges(3,4) * t557;
t719 = t553 * t559;
t716 = -qJ(5) - t556;
t705 = -t359 * qJ(4) - t495 * qJD(4);
t242 = t360 * pkin(3) - t705;
t482 = t495 * qJ(4);
t408 = pkin(3) * t496 + t482;
t713 = -t242 * t578 - t468 * t408;
t239 = rSges(7,1) * t342 + rSges(7,2) * t341 + rSges(7,3) * t494;
t712 = pkin(5) * t724 + t494 * t716 - t648 * t746 + t239;
t241 = t358 * pkin(3) + qJ(4) * t357 + qJD(4) * t497;
t382 = t469 * pkin(2) - pkin(9) * t468;
t363 = t555 * t382;
t711 = t555 * t241 + t363;
t683 = qJD(5) * t498;
t271 = t631 + t683;
t710 = -t241 - t271;
t383 = t471 * pkin(2) + t470 * pkin(9);
t709 = -t242 - t383;
t325 = -rSges(7,1) * t579 + rSges(7,2) * t475 + rSges(7,3) * t525;
t706 = pkin(5) * t721 + t716 * t525 + t659 * t746 + t325;
t362 = pkin(3) * t494 + qJ(4) * t493 + qJD(4) * t524;
t398 = rSges(5,1) * t648 - rSges(5,2) * t494 + rSges(5,3) * t493;
t704 = -t362 - t398;
t410 = pkin(4) * t648 + qJ(5) * t494 + qJD(5) * t525;
t703 = -t362 - t410;
t593 = rSges(5,1) * t576 - rSges(5,3) * t495;
t376 = -rSges(5,2) * t496 - t593;
t702 = -t376 - t408;
t377 = -rSges(5,1) * t578 - t498 * rSges(5,2) + t497 * rSges(5,3);
t409 = t498 * pkin(3) + qJ(4) * t497;
t701 = -t377 - t409;
t381 = t578 * t408;
t434 = -pkin(4) * t576 + t483;
t700 = -t434 * t578 - t381;
t399 = rSges(4,1) * t494 - rSges(4,2) * t493 + rSges(4,3) * t648;
t685 = qJD(2) * t553;
t516 = (pkin(2) * t733 + pkin(9) * t557) * t685;
t699 = -t399 - t516;
t400 = t409 * t648;
t698 = t435 * t648 + t400;
t472 = pkin(3) * t525 + qJ(4) * t524;
t697 = t408 * t659 - t472 * t576;
t474 = t529 * pkin(2) - pkin(9) * t578;
t460 = t555 * t474;
t696 = t555 * t409 + t460;
t695 = -t408 - t434;
t694 = -t409 - t435;
t530 = (pkin(2) * t557 - pkin(9) * t733) * t553;
t687 = qJD(1) * t558;
t650 = t553 * t687;
t503 = t530 * t650;
t693 = t472 * t650 + t503;
t445 = -rSges(5,1) * t659 - t525 * rSges(5,2) + t524 * rSges(5,3);
t692 = t445 + t472;
t446 = t525 * rSges(4,1) - t524 * rSges(4,2) - rSges(4,3) * t659;
t691 = -t446 - t530;
t473 = -pkin(2) * t577 - pkin(9) * t576;
t690 = t473 * t720 + t474 * t719;
t500 = -pkin(4) * t659 + t525 * qJ(5);
t689 = t472 + t500;
t686 = qJD(1) * t559;
t680 = m(6) / 0.2e1 + m(7) / 0.2e1;
t679 = -rSges(7,3) - pkin(3) + t556;
t676 = -t733 / 0.2e1;
t675 = t38 / 0.2e1 + t62 / 0.2e1;
t61 = t236 * t498 + t237 * t415 + t238 * t416 + t248 * t323 + t249 * t324 + t322 * t358;
t674 = t39 / 0.2e1 + t61 / 0.2e1;
t673 = t242 * t659 - t362 * t576 + t470 * t472;
t672 = t555 * t271 + t711;
t481 = t496 * qJD(5);
t272 = t481 + t745;
t671 = -t272 + t709;
t591 = -rSges(6,1) * t428 - rSges(6,2) * t427;
t294 = rSges(6,3) * t496 - t591;
t670 = -t294 + t695;
t295 = t430 * rSges(6,1) + t429 * rSges(6,2) + t498 * rSges(6,3);
t669 = -t295 + t694;
t176 = t304 * rSges(6,1) + t303 * rSges(6,2) + t358 * rSges(6,3);
t312 = rSges(6,1) * t437 + rSges(6,2) * t436 + rSges(6,3) * t494;
t668 = -t312 + t703;
t231 = t358 * rSges(4,1) - t357 * rSges(4,2) - t468 * rSges(4,3);
t649 = t553 * t686;
t666 = t382 * t719 + t383 * t720 + t473 * t649;
t361 = rSges(6,1) * t492 + rSges(6,2) * t491 + rSges(6,3) * t525;
t665 = t361 + t689;
t664 = -t516 + t704;
t663 = t555 * t435 + t696;
t662 = t500 * t650 + t693;
t661 = -t530 - t692;
t229 = -t468 * rSges(5,1) - t358 * rSges(5,2) + t357 * rSges(5,3);
t337 = t469 * rSges(3,1) + t468 * rSges(3,2) + rSges(3,3) * t649;
t379 = t498 * rSges(4,1) - t497 * rSges(4,2) - rSges(4,3) * t578;
t455 = t529 * rSges(3,1) + rSges(3,2) * t578 + rSges(3,3) * t720;
t653 = t733 * Icges(3,4);
t652 = t733 * t241;
t651 = t733 * t409;
t647 = t130 / 0.2e1 + t156 / 0.2e1;
t157 = t322 * t498 + t323 * t415 + t324 * t416;
t646 = t131 / 0.2e1 + t157 / 0.2e1;
t642 = -t558 * pkin(1) + pkin(8) * t719;
t641 = 2 * m(3);
t639 = 2 * m(4);
t637 = 2 * m(5);
t635 = 0.2e1 * m(6);
t633 = 0.2e1 * m(7);
t632 = t559 * t691;
t630 = -t272 * t578 - t468 * t434 + t713;
t629 = t703 - t712;
t628 = t695 - t708;
t627 = t694 - t707;
t626 = -t516 + t668;
t625 = t689 + t706;
t624 = -t530 - t665;
t623 = t408 * t720 + t409 * t719 + t690;
t622 = t434 * t659 - t500 * t576 + t697;
t24 = -t714 * t578 - t708 * t468 - (t710 - t715) * t576 + t627 * t470 + t630;
t592 = -t306 * rSges(6,1) - t305 * rSges(6,2);
t177 = t360 * rSges(6,3) - t592;
t37 = -t177 * t578 - t294 * t468 - (-t176 + t710) * t576 + t669 * t470 + t630;
t615 = m(6) * t37 + m(7) * t24;
t584 = t241 * t719 + t242 * t720 + t408 * t649 + t666;
t569 = t271 * t719 + t272 * t720 + t434 * t649 + t584;
t25 = (t715 * t559 + t714 * t558 + (t708 * t559 + (-t474 + t627) * t558) * qJD(1)) * t553 + t569;
t36 = (t176 * t559 + t177 * t558 + (t294 * t559 + (-t474 + t669) * t558) * qJD(1)) * t553 + t569;
t614 = m(6) * t36 + m(7) * t25;
t575 = -t271 * t733 - t652;
t30 = -t629 * t578 + t625 * t468 + (t684 * t707 - t715 * t733 + t575) * t553 + t698;
t54 = -t668 * t578 + t665 * t468 + (-t176 * t733 + t295 * t684 + t575) * t553 + t698;
t613 = m(6) * t54 + m(7) * t30;
t583 = t272 * t659 - t410 * t576 + t470 * t500 + t673;
t31 = -t712 * t576 + t706 * t470 + (t628 * t684 + t714 * t733) * t553 + t583;
t55 = -t576 * t312 + t470 * t361 + (t177 * t733 + t670 * t684) * t553 + t583;
t612 = m(6) * t55 + m(7) * t31;
t588 = -t516 + t629;
t40 = (t671 - t714) * t555 + (t559 * t588 + t687 * t706) * t553 + t662;
t68 = (-t177 + t671) * t555 + (t361 * t687 + t559 * t626) * t553 + t662;
t611 = m(6) * t68 + m(7) * t40;
t587 = -t530 - t625;
t581 = t587 * t559;
t41 = t715 * t555 + (qJD(1) * t581 + t558 * t588) * t553 + t672;
t585 = t624 * t559;
t69 = t555 * t176 + (qJD(1) * t585 + t558 * t626) * t553 + t672;
t610 = m(6) * t69 + m(7) * t41;
t605 = -pkin(1) * t687 + pkin(8) * t649;
t572 = t382 + t605;
t562 = t241 + t572;
t560 = t562 + t683;
t108 = t560 + t631 + t176;
t85 = t560 + t750;
t609 = m(6) * t108 + m(7) * t85;
t568 = -qJD(1) * t688 - t383;
t563 = t568 + t705;
t561 = -t481 + t563;
t109 = t360 * t734 + t561 + t592 - t745;
t86 = t360 * t679 - t470 * t547 + t561 + t590 - t677;
t608 = m(6) * t109 + m(7) * t86;
t582 = t434 * t720 + t435 * t719 + t623;
t119 = (t294 * t558 + t295 * t559) * t553 + t582;
t91 = (t558 * t708 + t559 * t707) * t553 + t582;
t607 = m(6) * t119 + m(7) * t91;
t120 = -t294 * t578 - t576 * t669 + t700;
t92 = -t576 * t627 - t578 * t708 + t700;
t606 = m(6) * t120 + m(7) * t92;
t604 = t661 * t559;
t110 = t555 * t707 + t587 * t720 + t663;
t136 = t295 * t555 + t624 * t720 + t663;
t603 = m(6) * t136 + m(7) * t110;
t111 = (-t473 + t628) * t555 + t553 * t581;
t137 = (-t473 + t670) * t555 + t553 * t585;
t602 = m(6) * t137 + m(7) * t111;
t112 = -t576 * t706 + t659 * t708 + t622;
t150 = t294 * t659 - t361 * t576 + t622;
t601 = m(6) * t150 + m(7) * t112;
t574 = -t435 * t733 - t651;
t113 = (-t707 * t733 + t574) * t553 + t625 * t578;
t151 = (-t295 * t733 + t574) * t553 + t665 * t578;
t600 = m(6) * t151 + m(7) * t113;
t580 = -t473 + t642;
t573 = -t482 + t580;
t168 = t496 * t679 + t547 * t576 + t573 + t589 - t678;
t194 = t496 * t734 - t434 + t573 + t591;
t599 = m(6) * t194 + m(7) * t168;
t586 = t474 + t688;
t570 = t409 + t586;
t169 = t570 + t749;
t195 = t435 + t570 + t295;
t598 = m(6) * t195 + m(7) * t169;
t597 = -t198 / 0.2e1 - t200 / 0.2e1 - t132 / 0.2e1;
t596 = t199 / 0.2e1 + t201 / 0.2e1 + t133 / 0.2e1;
t595 = -t471 * rSges(3,1) + t470 * rSges(3,2);
t594 = -t470 * rSges(5,1) - t359 * rSges(5,3);
t232 = t360 * rSges(4,1) - t359 * rSges(4,2) + t470 * rSges(4,3);
t378 = rSges(4,1) * t496 - rSges(4,2) * t495 - rSges(4,3) * t576;
t454 = -rSges(3,1) * t577 + rSges(3,2) * t576 - rSges(3,3) * t719;
t507 = Icges(3,6) * t555 + (Icges(3,2) * t733 + t727) * t553;
t508 = Icges(3,5) * t555 + (Icges(3,1) * t557 + t653) * t553;
t512 = (Icges(3,5) * t733 - Icges(3,6) * t557) * t685;
t513 = (-Icges(3,2) * t557 + t653) * t685;
t514 = (Icges(3,1) * t733 - t727) * t685;
t571 = t553 * t557 * t514 - t507 * t648 + t508 * t617 + t555 * t512 + t513 * t659;
t116 = t359 * t438 - t360 * t439 + t392 * t495 - t393 * t496 - t394 * t576 + t440 * t470;
t118 = -t359 * t442 + t360 * t443 - t395 * t576 - t396 * t495 + t397 * t496 + t441 * t470;
t82 = t305 * t355 + t306 * t356 + t309 * t496 + t310 * t427 + t311 * t428 + t354 * t360;
t567 = t89 / 0.2e1 + t87 / 0.2e1 + t48 / 0.2e1 + t118 / 0.2e1 + t116 / 0.2e1 + t82 / 0.2e1 + t675;
t115 = t357 * t438 - t358 * t439 + t392 * t497 - t393 * t498 - t394 * t578 - t440 * t468;
t117 = -t357 * t442 + t358 * t443 - t395 * t578 - t396 * t497 + t397 * t498 - t441 * t468;
t81 = t303 * t355 + t304 * t356 + t309 * t498 + t310 * t429 + t311 * t430 + t354 * t358;
t566 = t90 / 0.2e1 + t88 / 0.2e1 + t49 / 0.2e1 + t117 / 0.2e1 + t115 / 0.2e1 + t81 / 0.2e1 + t674;
t160 = t354 * t496 + t355 * t427 + t356 * t428;
t213 = t438 * t495 - t439 * t496 - t440 * t576;
t215 = -t441 * t576 - t442 * t495 + t443 * t496;
t565 = t213 / 0.2e1 + t215 / 0.2e1 + t160 / 0.2e1 - t597 + t647;
t161 = t354 * t498 + t355 * t429 + t356 * t430;
t214 = t438 * t497 - t439 * t498 - t440 * t578;
t216 = -t441 * t578 - t442 * t497 + t443 * t498;
t564 = -t216 / 0.2e1 - t214 / 0.2e1 - t161 / 0.2e1 - t596 - t646;
t515 = (rSges(3,1) * t733 - rSges(3,2) * t557) * t685;
t511 = t555 * rSges(3,3) + (rSges(3,1) * t557 + rSges(3,2) * t733) * t553;
t506 = Icges(3,3) * t555 + (Icges(3,5) * t557 + Icges(3,6) * t733) * t553;
t452 = Icges(3,1) * t529 + Icges(3,4) * t578 + Icges(3,5) * t720;
t451 = -Icges(3,1) * t577 + Icges(3,4) * t576 - Icges(3,5) * t719;
t450 = Icges(3,4) * t529 + Icges(3,2) * t578 + Icges(3,6) * t720;
t449 = -Icges(3,4) * t577 + Icges(3,2) * t576 - Icges(3,6) * t719;
t448 = Icges(3,5) * t529 + Icges(3,6) * t578 + Icges(3,3) * t720;
t447 = -Icges(3,5) * t577 + Icges(3,6) * t576 - Icges(3,3) * t719;
t423 = t455 + t688;
t422 = -t454 + t642;
t390 = -t555 * t454 - t511 * t719;
t389 = t455 * t555 - t511 * t720;
t338 = rSges(3,3) * t650 - t595;
t335 = Icges(3,1) * t471 - Icges(3,4) * t470 + Icges(3,5) * t650;
t334 = Icges(3,1) * t469 + Icges(3,4) * t468 + Icges(3,5) * t649;
t333 = Icges(3,4) * t471 - Icges(3,2) * t470 + Icges(3,6) * t650;
t332 = Icges(3,4) * t469 + Icges(3,2) * t468 + Icges(3,6) * t649;
t331 = Icges(3,5) * t471 - Icges(3,6) * t470 + Icges(3,3) * t650;
t330 = Icges(3,5) * t469 + Icges(3,6) * t468 + Icges(3,3) * t649;
t319 = (-t550 + (-rSges(3,3) - pkin(8)) * t720) * qJD(1) + t595;
t318 = t605 + t337;
t317 = t506 * t720 + t507 * t578 + t508 * t529;
t316 = -t506 * t719 + t507 * t576 - t508 * t577;
t308 = t586 + t379;
t307 = -t378 + t580;
t298 = t571 * t555;
t286 = t555 * t337 + (-t511 * t686 - t515 * t558) * t553;
t285 = -t555 * t338 + (t511 * t687 - t515 * t559) * t553;
t284 = -t379 * t659 + t446 * t578;
t283 = t378 * t659 - t446 * t576;
t280 = t555 * t448 + (t450 * t733 + t452 * t557) * t553;
t279 = t555 * t447 + (t449 * t733 + t451 * t557) * t553;
t263 = t448 * t720 + t450 * t578 + t452 * t529;
t262 = t447 * t720 + t449 * t578 + t451 * t529;
t261 = -t448 * t719 + t450 * t576 - t452 * t577;
t260 = -t447 * t719 + t449 * t576 - t451 * t577;
t259 = -t441 * t659 - t524 * t442 + t525 * t443;
t258 = t524 * t438 - t525 * t439 - t440 * t659;
t255 = t259 * t648;
t254 = t258 * t648;
t253 = t570 + t377;
t252 = t496 * t735 + t573 + t593;
t245 = -t378 * t578 + t379 * t576;
t244 = (-t378 - t473) * t555 + t553 * t632;
t243 = t379 * t555 + t691 * t720 + t460;
t230 = -t360 * rSges(5,2) - t594;
t208 = -t470 * t507 + t471 * t508 + t576 * t513 - t577 * t514 + (t506 * t687 - t512 * t559) * t553;
t207 = t468 * t507 + t469 * t508 + t578 * t513 + t529 * t514 + (t506 * t686 + t512 * t558) * t553;
t206 = (t378 * t558 + t379 * t559) * t553 + t690;
t205 = (-t377 * t733 - t651) * t553 + t692 * t578;
t204 = t376 * t659 - t445 * t576 + t697;
t197 = t282 * t525 - t325 * t498;
t196 = -t281 * t525 + t325 * t496;
t193 = -t232 + t568;
t192 = t572 + t231;
t191 = (-t473 + t702) * t555 + t553 * t604;
t190 = t377 * t555 + t661 * t720 + t696;
t189 = -t367 * t578 - t371 * t497 + t375 * t498;
t188 = -t366 * t578 - t370 * t497 + t374 * t498;
t187 = -t367 * t576 - t371 * t495 + t375 * t496;
t186 = -t366 * t576 - t370 * t495 + t374 * t496;
t185 = t365 * t497 - t369 * t498 - t373 * t578;
t184 = t364 * t497 - t368 * t498 - t372 * t578;
t183 = t365 * t495 - t369 * t496 - t373 * t576;
t182 = t364 * t495 - t368 * t496 - t372 * t576;
t181 = t354 * t525 + t355 * t491 + t356 * t492;
t180 = t181 * t648;
t167 = -t376 * t578 - t576 * t701 - t381;
t166 = t281 * t498 - t282 * t496;
t164 = t165 * t648;
t163 = t555 * t330 + (t733 * t332 + t334 * t557 + (-t450 * t557 + t452 * t733) * qJD(2)) * t553;
t162 = t555 * t331 + (t733 * t333 + t335 * t557 + (-t449 * t557 + t451 * t733) * qJD(2)) * t553;
t158 = (t376 * t558 + t377 * t559) * t553 + t623;
t141 = t555 * t231 + t363 + (qJD(1) * t632 + t558 * t699) * t553;
t140 = t503 + (-t232 - t383) * t555 + (t446 * t687 + t559 * t699) * t553;
t139 = -t576 * t399 + t470 * t446 + (t232 * t733 - t378 * t684) * t553;
t138 = t578 * t399 + t468 * t446 + (-t231 * t733 + t379 * t684) * t553;
t135 = t360 * t735 + t563 + t594;
t134 = t562 + t229;
t128 = t289 * t498 + t291 * t429 + t293 * t430;
t127 = t288 * t498 + t290 * t429 + t292 * t430;
t126 = t289 * t496 + t291 * t427 + t293 * t428;
t125 = t288 * t496 + t290 * t427 + t292 * t428;
t124 = t274 * t498 + t276 * t415 + t278 * t416;
t123 = t273 * t498 + t275 * t415 + t277 * t416;
t114 = t231 * t576 - t232 * t578 - t378 * t468 - t379 * t470;
t107 = t555 * t229 + (qJD(1) * t604 + t558 * t664) * t553 + t711;
t106 = (-t230 + t709) * t555 + (t445 * t687 + t559 * t664) * t553 + t693;
t103 = (t231 * t559 + t232 * t558 + (t378 * t559 + (-t379 - t474) * t558) * qJD(1)) * t553 + t666;
t102 = t216 * t555 + (-t188 * t559 + t189 * t558) * t553;
t101 = t215 * t555 + (-t186 * t559 + t187 * t558) * t553;
t100 = t214 * t555 + (-t184 * t559 + t185 * t558) * t553;
t99 = t213 * t555 + (-t182 * t559 + t183 * t558) * t553;
t98 = -t188 * t576 - t189 * t578 - t216 * t659;
t97 = -t186 * t576 - t187 * t578 - t215 * t659;
t96 = -t184 * t576 - t185 * t578 - t214 * t659;
t95 = -t182 * t576 - t183 * t578 - t213 * t659;
t94 = -t576 * t398 + t470 * t445 + (t230 * t733 + t684 * t702) * t553 + t673;
t93 = t400 - t704 * t578 + t692 * t468 + (-t229 * t733 + t377 * t684 - t652) * t553;
t84 = -t149 * t525 + t239 * t496 - t281 * t494 + t325 * t360;
t83 = t148 * t525 - t239 * t498 + t282 * t494 - t325 * t358;
t80 = -t219 * t576 - t223 * t495 + t227 * t496 - t359 * t371 + t360 * t375 + t367 * t470;
t79 = -t220 * t576 - t224 * t495 + t228 * t496 - t359 * t370 + t360 * t374 + t366 * t470;
t78 = -t219 * t578 - t223 * t497 + t227 * t498 - t357 * t371 + t358 * t375 - t367 * t468;
t77 = -t220 * t578 - t224 * t497 + t228 * t498 - t357 * t370 + t358 * t374 - t366 * t468;
t76 = t217 * t495 - t221 * t496 - t225 * t576 + t359 * t365 - t360 * t369 + t373 * t470;
t75 = t218 * t495 - t222 * t496 - t226 * t576 + t359 * t364 - t360 * t368 + t372 * t470;
t74 = t217 * t497 - t221 * t498 - t225 * t578 + t357 * t365 - t358 * t369 - t373 * t468;
t73 = t218 * t497 - t222 * t498 - t226 * t578 + t357 * t364 - t358 * t368 - t372 * t468;
t67 = -t230 * t578 - t376 * t468 - (-t229 - t241) * t576 + t701 * t470 + t713;
t66 = (t229 * t559 + t230 * t558 + (t376 * t559 + (-t474 + t701) * t558) * qJD(1)) * t553 + t584;
t65 = t165 * t555 + (-t130 * t559 + t131 * t558) * t553;
t64 = -t148 * t496 + t149 * t498 + t281 * t358 - t282 * t360;
t63 = -t130 * t576 - t131 * t578 - t165 * t659;
t60 = t130 * t496 + t131 * t498 + t165 * t525;
t59 = t161 * t555 + (-t127 * t559 + t128 * t558) * t553;
t58 = t160 * t555 + (-t125 * t559 + t126 * t558) * t553;
t57 = -t127 * t576 - t128 * t578 - t161 * t659;
t56 = -t125 * t576 - t126 * t578 - t160 * t659;
t53 = t157 * t555 + (-t123 * t559 + t124 * t558) * t553;
t52 = t156 * t555 + (-t121 * t559 + t122 * t558) * t553;
t51 = -t123 * t576 - t124 * t578 - t157 * t659;
t50 = -t121 * t576 - t122 * t578 - t156 * t659;
t47 = t123 * t496 + t124 * t498 + t157 * t525;
t46 = t121 * t496 + t122 * t498 + t156 * t525;
t45 = t170 * t496 + t172 * t427 + t174 * t428 + t289 * t360 + t291 * t305 + t293 * t306;
t44 = t171 * t496 + t173 * t427 + t175 * t428 + t288 * t360 + t290 * t305 + t292 * t306;
t43 = t170 * t498 + t172 * t429 + t174 * t430 + t289 * t358 + t291 * t303 + t293 * t304;
t42 = t171 * t498 + t173 * t429 + t175 * t430 + t288 * t358 + t290 * t303 + t292 * t304;
t33 = t142 * t498 + t144 * t415 + t146 * t416 + t248 * t276 + t249 * t278 + t274 * t358;
t32 = t143 * t498 + t145 * t415 + t147 * t416 + t248 * t275 + t249 * t277 + t273 * t358;
t27 = -t155 * t659 + t200 * t470 - t201 * t468 - t576 * t89 - t578 * t90 + t255;
t26 = -t154 * t659 + t198 * t470 - t199 * t468 - t576 * t87 - t578 * t88 + t254;
t23 = t118 * t555 + (t558 * t80 - t559 * t79 + (t186 * t558 + t187 * t559) * qJD(1)) * t553;
t22 = t117 * t555 + (t558 * t78 - t559 * t77 + (t188 * t558 + t189 * t559) * qJD(1)) * t553;
t21 = t116 * t555 + (t558 * t76 - t559 * t75 + (t182 * t558 + t183 * t559) * qJD(1)) * t553;
t20 = t115 * t555 + (t558 * t74 - t559 * t73 + (t184 * t558 + t185 * t559) * qJD(1)) * t553;
t19 = t186 * t470 - t187 * t468 - t79 * t576 - t80 * t578 + (-t118 * t733 + t215 * t684) * t553;
t18 = t188 * t470 - t189 * t468 - t77 * t576 - t78 * t578 + (-t117 * t733 + t216 * t684) * t553;
t17 = t182 * t470 - t183 * t468 - t75 * t576 - t76 * t578 + (-t116 * t733 + t213 * t684) * t553;
t16 = t184 * t470 - t185 * t468 - t73 * t576 - t74 * t578 + (-t115 * t733 + t214 * t684) * t553;
t14 = -t105 * t659 + t132 * t470 - t133 * t468 - t48 * t576 - t49 * t578 + t180;
t13 = t82 * t555 + (-t44 * t559 + t45 * t558 + (t125 * t558 + t126 * t559) * qJD(1)) * t553;
t12 = t81 * t555 + (-t42 * t559 + t43 * t558 + (t127 * t558 + t128 * t559) * qJD(1)) * t553;
t10 = t125 * t470 - t126 * t468 - t44 * t576 - t45 * t578 + (t160 * t684 - t733 * t82) * t553;
t9 = t127 * t470 - t128 * t468 - t42 * t576 - t43 * t578 + (t161 * t684 - t733 * t81) * t553;
t8 = t130 * t470 - t131 * t468 - t38 * t576 - t39 * t578 - t659 * t72 + t164;
t7 = t130 * t360 + t131 * t358 + t38 * t496 + t39 * t498 + t728;
t6 = t62 * t555 + (-t34 * t559 + t35 * t558 + (t121 * t558 + t122 * t559) * qJD(1)) * t553;
t5 = t61 * t555 + (-t32 * t559 + t33 * t558 + (t123 * t558 + t124 * t559) * qJD(1)) * t553;
t4 = t121 * t470 - t122 * t468 - t34 * t576 - t35 * t578 + (t156 * t684 - t62 * t733) * t553;
t3 = t123 * t470 - t124 * t468 - t32 * t576 - t33 * t578 + (t157 * t684 - t61 * t733) * t553;
t1 = t123 * t360 + t124 * t358 + t157 * t494 + t32 * t496 + t33 * t498 + t525 * t61;
t15 = [(t108 * t195 + t109 * t194) * t635 + (t134 * t253 + t135 * t252) * t637 + (t168 * t86 + t169 * t85) * t633 + (t192 * t308 + t193 * t307) * t639 + (t318 * t423 + t319 * t422) * t641 + t571 - t743; (t140 * t307 + t141 * t308 + t192 * t243 + t193 * t244) * m(4) + (t106 * t252 + t107 * t253 + t134 * t190 + t135 * t191) * m(5) + (t108 * t136 + t109 * t137 + t194 * t68 + t195 * t69) * m(6) + (t110 * t85 + t111 * t86 + t168 * t40 + t169 * t41) * m(7) + t71 + t298 + m(3) * (t285 * t422 + t286 * t423 + t318 * t389 + t319 * t390) + ((-t208 / 0.2e1 - t162 / 0.2e1 - t567) * t559 + (t207 / 0.2e1 + t163 / 0.2e1 + t566) * t558 + ((t317 / 0.2e1 + t280 / 0.2e1 - t564) * t559 + (t316 / 0.2e1 + t279 / 0.2e1 + t565) * t558) * qJD(1)) * t553 + t748; (t110 * t41 + t111 * t40 + t25 * t91) * t633 + (t119 * t36 + t136 * t69 + t137 * t68) * t635 + (t106 * t191 + t107 * t190 + t158 * t66) * t637 + (t103 * t206 + t140 * t244 + t141 * t243) * t639 + (t390 * t285 + t389 * t286 + (t454 * t558 + t455 * t559) * (t337 * t559 + t338 * t558 + (t454 * t559 - t455 * t558) * qJD(1)) * t553 ^ 2) * t641 + (((t332 * t578 + t529 * t334 + t468 * t450 + t469 * t452) * t558 + t263 * t686 - (t333 * t578 + t529 * t335 + t468 * t449 + t469 * t451) * t559 + t262 * t687 + ((t330 * t558 + t448 * t686) * t558 - (t331 * t558 + t447 * t686) * t559) * t553) * t553 + t22 + t20 + t12 + t5) * t720 + (-t23 - t21 - t13 - t6 - ((t332 * t576 - t334 * t577 - t470 * t450 + t471 * t452) * t558 + t261 * t686 - (t333 * t576 - t335 * t577 - t470 * t449 + t471 * t451) * t559 + t260 * t687 + ((-t330 * t559 + t448 * t687) * t558 - (-t331 * t559 + t447 * t687) * t559) * t553) * t553) * t719 + (t101 + t99 + t58 + t52 + (-t260 * t559 + t261 * t558) * t553) * t650 + (t102 + t100 + t59 + t53 + (-t262 * t559 + t263 * t558) * t553) * t649 + (-t208 * t719 + t317 * t649 + t316 * t650 + t298 + (-t162 * t559 + t163 * t558 + (t279 * t558 + t280 * t559) * qJD(1)) * t553 + t207 * t720 + t744) * t555; (t138 * t308 + t139 * t307 + t192 * t284 + t193 * t283) * m(4) + (t134 * t205 + t135 * t204 + t252 * t94 + t253 * t93) * m(5) + (t108 * t151 + t109 * t150 + t194 * t55 + t195 * t54) * m(6) + (t112 * t86 + t113 * t85 + t168 * t31 + t169 * t30) * m(7) + t164 + t180 + t254 + t255 - t566 * t578 - t567 * t576 + t565 * t470 + t564 * t468 + t743 * t659; (t103 * t245 + t114 * t206 + t138 * t243 + t139 * t244 + t140 * t283 + t141 * t284) * m(4) + (t106 * t204 + t107 * t205 + t158 * t67 + t167 * t66 + t190 * t93 + t191 * t94) * m(5) + (t119 * t37 + t120 * t36 + t136 * t54 + t137 * t55 + t150 * t68 + t151 * t69) * m(6) + (t110 * t30 + t111 * t31 + t112 * t40 + t113 * t41 + t24 * t91 + t25 * t92) * m(7) + (t27 / 0.2e1 + t26 / 0.2e1 + t14 / 0.2e1 + t8 / 0.2e1) * t555 - (t22 / 0.2e1 + t20 / 0.2e1 + t12 / 0.2e1 + t5 / 0.2e1) * t578 - (t23 / 0.2e1 + t21 / 0.2e1 + t13 / 0.2e1 + t6 / 0.2e1) * t576 + (t101 / 0.2e1 + t99 / 0.2e1 + t58 / 0.2e1 + t52 / 0.2e1) * t470 + (-t100 / 0.2e1 - t102 / 0.2e1 - t59 / 0.2e1 - t53 / 0.2e1) * t468 + ((-t4 / 0.2e1 - t10 / 0.2e1 - t17 / 0.2e1 - t19 / 0.2e1) * t559 + (t3 / 0.2e1 + t9 / 0.2e1 + t16 / 0.2e1 + t18 / 0.2e1) * t558 + (t65 / 0.2e1 + t597 * t719 + t596 * t720 + (t258 / 0.2e1 + t259 / 0.2e1 + t181 / 0.2e1) * t555) * t684 + ((t96 / 0.2e1 + t98 / 0.2e1 + t51 / 0.2e1 + t57 / 0.2e1) * t559 + (t95 / 0.2e1 + t97 / 0.2e1 + t50 / 0.2e1 + t56 / 0.2e1) * t558) * qJD(1) + t744 * t676) * t553; (t114 * t245 + t138 * t284 + t139 * t283) * t639 + (t167 * t67 + t204 * t94 + t205 * t93) * t637 + (t120 * t37 + t150 * t55 + t151 * t54) * t635 + (t112 * t31 + t113 * t30 + t24 * t92) * t633 + (-t14 - t26 - t27 - t8) * t659 - (t18 + t16 + t9 + t3) * t578 - (t19 + t17 + t10 + t4) * t576 + (t97 + t95 + t56 + t50) * t470 + (-t98 - t96 - t57 - t51) * t468 + (t63 - t751 * t578 - t752 * t576 + (-t181 - t258 - t259) * t659) * t648; (m(5) * t135 + t608) * t497 + (m(5) * t134 + t609) * t495 + (m(5) * t253 + t598) * t359 + (m(5) * t252 + t599) * t357; (m(5) * t66 + t614) * t524 + (m(5) * t106 + t611) * t497 + (m(5) * t107 + t610) * t495 + (m(5) * t158 + t607) * t493 + (m(5) * t190 + t603) * t359 + (m(5) * t191 + t602) * t357; (m(5) * t67 + t615) * t524 + (m(5) * t94 + t612) * t497 + (m(5) * t93 + t613) * t495 + (m(5) * t167 + t606) * t493 + (m(5) * t205 + t600) * t359 + (m(5) * t204 + t601) * t357; 0.4e1 * (m(5) / 0.2e1 + t680) * (t357 * t497 + t359 * t495 + t493 * t524); t358 * t599 + t360 * t598 + t496 * t609 + t498 * t608; t358 * t602 + t360 * t603 + t494 * t607 + t496 * t610 + t498 * t611 + t525 * t614; t358 * t601 + t360 * t600 + t494 * t606 + t496 * t613 + t498 * t612 + t525 * t615; 0.2e1 * t680 * (t357 * t498 + t358 * t497 + t359 * t496 + t360 * t495 + t493 * t525 + t494 * t524); 0.4e1 * t680 * (t358 * t498 + t360 * t496 + t494 * t525); (t168 * t84 + t169 * t83 + t196 * t86 + t197 * t85) * m(7) + t674 * t498 + t675 * t496 + t647 * t360 + t646 * t358 + t728; (t110 * t83 + t111 * t84 + t166 * t25 + t196 * t40 + t197 * t41 + t64 * t91) * m(7) + t53 * t742 + t5 * t738 + t555 * t7 / 0.2e1 + t52 * t741 + t6 * t739 + t65 * t740 + t11 * t737 + (t1 * t736 + t559 * t747 + (t559 * t47 / 0.2e1 + t46 * t736) * qJD(1)) * t553; t63 * t740 + t8 * t737 + t51 * t742 + t3 * t738 - t468 * t47 / 0.2e1 - t578 * t1 / 0.2e1 + t470 * t46 / 0.2e1 + t576 * t747 + (t112 * t84 + t113 * t83 + t166 * t24 + t196 * t31 + t197 * t30 + t64 * t92) * m(7) + t50 * t741 + t4 * t739 + (t60 * t684 / 0.2e1 + t7 * t676) * t553; (t166 * t493 + t196 * t357 + t197 * t359 + t495 * t83 + t497 * t84 + t524 * t64) * m(7); (t166 * t494 + t196 * t358 + t197 * t360 + t496 * t83 + t498 * t84 + t525 * t64) * m(7); t358 * t47 + t498 * t1 + t360 * t46 + t496 * t2 + t494 * t60 + t525 * t7 + (t166 * t64 + t196 * t84 + t197 * t83) * t633;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
