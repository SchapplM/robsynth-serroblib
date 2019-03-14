% Calculate time derivative of joint inertia matrix for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR12_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:18
% EndTime: 2019-03-09 14:36:19
% DurationCPUTime: 33.89s
% Computational Cost: add. (119747->1451), mult. (219843->1913), div. (0->0), fcn. (243419->12), ass. (0->605)
t618 = cos(pkin(6));
t622 = sin(qJ(1));
t625 = cos(qJ(2));
t774 = t622 * t625;
t621 = sin(qJ(2));
t626 = cos(qJ(1));
t775 = t621 * t626;
t587 = t618 * t774 + t775;
t776 = t621 * t622;
t723 = t618 * t776;
t773 = t625 * t626;
t588 = -t723 + t773;
t617 = sin(pkin(6));
t781 = t617 * t622;
t480 = Icges(4,1) * t781 - Icges(4,4) * t588 + Icges(4,5) * t587;
t483 = Icges(3,5) * t588 - Icges(3,6) * t587 + Icges(3,3) * t781;
t833 = t480 + t483;
t722 = t618 * t773;
t585 = -t722 + t776;
t586 = t618 * t775 + t774;
t779 = t617 * t626;
t481 = -Icges(4,1) * t779 - Icges(4,4) * t586 + Icges(4,5) * t585;
t482 = Icges(3,5) * t586 - Icges(3,6) * t585 - Icges(3,3) * t779;
t832 = t481 + t482;
t477 = -Icges(4,5) * t779 - Icges(4,6) * t586 + Icges(4,3) * t585;
t484 = Icges(3,4) * t586 - Icges(3,2) * t585 - Icges(3,6) * t779;
t831 = -t484 + t477;
t476 = Icges(4,5) * t781 - Icges(4,6) * t588 + Icges(4,3) * t587;
t485 = Icges(3,4) * t588 - Icges(3,2) * t587 + Icges(3,6) * t781;
t830 = t485 - t476;
t479 = -Icges(4,4) * t779 - Icges(4,2) * t586 + Icges(4,6) * t585;
t486 = Icges(3,1) * t586 - Icges(3,4) * t585 - Icges(3,5) * t779;
t829 = t486 - t479;
t478 = Icges(4,4) * t781 - Icges(4,2) * t588 + Icges(4,6) * t587;
t487 = Icges(3,1) * t588 - Icges(3,4) * t587 + Icges(3,5) * t781;
t828 = -t487 + t478;
t796 = Icges(3,4) * t621;
t553 = Icges(3,6) * t618 + (Icges(3,2) * t625 + t796) * t617;
t795 = Icges(3,4) * t625;
t554 = Icges(3,5) * t618 + (Icges(3,1) * t621 + t795) * t617;
t794 = Icges(4,6) * t621;
t555 = Icges(4,5) * t618 + (-Icges(4,3) * t625 - t794) * t617;
t735 = qJD(2) * t617;
t567 = (Icges(3,5) * t625 - Icges(3,6) * t621) * t735;
t568 = (-Icges(4,4) * t625 + Icges(4,5) * t621) * t735;
t569 = (-Icges(3,2) * t621 + t795) * t735;
t570 = (Icges(3,1) * t625 - t796) * t735;
t734 = qJD(2) * t621;
t704 = t617 * t734;
t733 = qJD(2) * t625;
t705 = t617 * t733;
t780 = t617 * t625;
t782 = t617 * t621;
t857 = t554 * t705 + t569 * t780 + t570 * t782 + (-t553 + t555) * t704 + (t567 + t568) * t618;
t772 = qJ(4) + qJ(5);
t614 = sin(t772);
t694 = cos(t772);
t664 = t617 * t694;
t526 = t587 * t614 + t622 * t664;
t619 = sin(qJ(6));
t623 = cos(qJ(6));
t453 = -t526 * t619 + t588 * t623;
t454 = t526 * t623 + t588 * t619;
t525 = -t587 * t694 + t614 * t781;
t309 = Icges(7,5) * t454 + Icges(7,6) * t453 + Icges(7,3) * t525;
t311 = Icges(7,4) * t454 + Icges(7,2) * t453 + Icges(7,6) * t525;
t313 = Icges(7,1) * t454 + Icges(7,4) * t453 + Icges(7,5) * t525;
t528 = t585 * t614 - t626 * t664;
t455 = -t528 * t619 + t586 * t623;
t456 = t528 * t623 + t586 * t619;
t527 = t585 * t694 + t614 * t779;
t150 = -t309 * t527 + t311 * t455 + t313 * t456;
t310 = Icges(7,5) * t456 + Icges(7,6) * t455 - Icges(7,3) * t527;
t312 = Icges(7,4) * t456 + Icges(7,2) * t455 - Icges(7,6) * t527;
t314 = Icges(7,1) * t456 + Icges(7,4) * t455 - Icges(7,5) * t527;
t151 = -t310 * t527 + t312 * t455 + t314 * t456;
t562 = -t614 * t780 + t618 * t694;
t522 = -t562 * t619 + t623 * t782;
t523 = t562 * t623 + t619 * t782;
t561 = t618 * t614 + t625 * t664;
t366 = Icges(7,5) * t523 + Icges(7,6) * t522 + Icges(7,3) * t561;
t367 = Icges(7,4) * t523 + Icges(7,2) * t522 + Icges(7,6) * t561;
t368 = Icges(7,1) * t523 + Icges(7,4) * t522 + Icges(7,5) * t561;
t188 = -t366 * t527 + t367 * t455 + t368 * t456;
t678 = qJD(2) * t618 + qJD(1);
t512 = -qJD(1) * t722 - t626 * t733 + t678 * t776;
t730 = qJD(4) + qJD(5);
t647 = t730 * t694;
t653 = qJD(1) * t664;
t681 = t617 * t730;
t363 = t626 * t653 + t587 * t647 + (-t622 * t681 - t512) * t614;
t513 = qJD(1) * t586 + qJD(2) * t587;
t274 = -qJD(6) * t454 - t363 * t619 - t513 * t623;
t275 = qJD(6) * t453 + t363 * t623 - t513 * t619;
t642 = t617 * t647;
t736 = qJD(1) * t626;
t707 = t617 * t736;
t362 = t622 * t642 + t512 * t694 + (t587 * t730 + t707) * t614;
t171 = Icges(7,5) * t275 + Icges(7,6) * t274 + Icges(7,3) * t362;
t173 = Icges(7,4) * t275 + Icges(7,2) * t274 + Icges(7,6) * t362;
t175 = Icges(7,1) * t275 + Icges(7,4) * t274 + Icges(7,5) * t362;
t514 = qJD(1) * t587 + qJD(2) * t586;
t361 = t622 * t653 + t585 * t647 + (t626 * t681 + t514) * t614;
t515 = -qJD(1) * t723 - t622 * t734 + t678 * t773;
t272 = -qJD(6) * t456 - t361 * t619 + t515 * t623;
t273 = qJD(6) * t455 + t361 * t623 + t515 * t619;
t737 = qJD(1) * t622;
t708 = t617 * t737;
t360 = -t626 * t642 - t514 * t694 + (t585 * t730 + t708) * t614;
t48 = -t171 * t527 + t173 * t455 + t175 * t456 + t272 * t311 + t273 * t313 + t309 * t360;
t170 = Icges(7,5) * t273 + Icges(7,6) * t272 + Icges(7,3) * t360;
t172 = Icges(7,4) * t273 + Icges(7,2) * t272 + Icges(7,6) * t360;
t174 = Icges(7,1) * t273 + Icges(7,4) * t272 + Icges(7,5) * t360;
t49 = -t170 * t527 + t172 * t455 + t174 * t456 + t272 * t312 + t273 * t314 + t310 * t360;
t500 = -t625 * t642 + (-t618 * t730 + t704) * t614;
t378 = -qJD(6) * t523 - t500 * t619 + t623 * t705;
t379 = qJD(6) * t522 + t500 * t623 + t619 * t705;
t501 = -t614 * t625 * t681 + t618 * t647 - t664 * t734;
t251 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t501;
t252 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t501;
t253 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t501;
t80 = -t251 * t527 + t252 * t455 + t253 * t456 + t272 * t367 + t273 * t368 + t360 * t366;
t11 = -t150 * t513 + t151 * t515 + t48 * t588 + t49 * t586 + (t188 * t733 + t621 * t80) * t617;
t369 = Icges(6,5) * t500 - Icges(6,6) * t501 + Icges(6,3) * t705;
t370 = Icges(6,4) * t500 - Icges(6,2) * t501 + Icges(6,6) * t705;
t371 = Icges(6,1) * t500 - Icges(6,4) * t501 + Icges(6,5) * t705;
t461 = Icges(6,5) * t562 - Icges(6,6) * t561 + Icges(6,3) * t782;
t462 = Icges(6,4) * t562 - Icges(6,2) * t561 + Icges(6,6) * t782;
t463 = Icges(6,1) * t562 - Icges(6,4) * t561 + Icges(6,5) * t782;
t139 = -t360 * t462 + t361 * t463 + t369 * t586 + t370 * t527 + t371 * t528 + t461 * t515;
t383 = Icges(6,5) * t526 - Icges(6,6) * t525 + Icges(6,3) * t588;
t385 = Icges(6,4) * t526 - Icges(6,2) * t525 + Icges(6,6) * t588;
t387 = Icges(6,1) * t526 - Icges(6,4) * t525 + Icges(6,5) * t588;
t208 = t383 * t586 + t385 * t527 + t387 * t528;
t384 = Icges(6,5) * t528 + Icges(6,6) * t527 + Icges(6,3) * t586;
t386 = Icges(6,4) * t528 + Icges(6,2) * t527 + Icges(6,6) * t586;
t388 = Icges(6,1) * t528 + Icges(6,4) * t527 + Icges(6,5) * t586;
t209 = t384 * t586 + t386 * t527 + t388 * t528;
t240 = t461 * t586 + t462 * t527 + t463 * t528;
t244 = Icges(6,5) * t363 - Icges(6,6) * t362 - Icges(6,3) * t513;
t246 = Icges(6,4) * t363 - Icges(6,2) * t362 - Icges(6,6) * t513;
t248 = Icges(6,1) * t363 - Icges(6,4) * t362 - Icges(6,5) * t513;
t94 = t244 * t586 + t246 * t527 + t248 * t528 - t360 * t385 + t361 * t387 + t383 * t515;
t243 = Icges(6,5) * t361 - Icges(6,6) * t360 + Icges(6,3) * t515;
t245 = Icges(6,4) * t361 - Icges(6,2) * t360 + Icges(6,6) * t515;
t247 = Icges(6,1) * t361 - Icges(6,4) * t360 + Icges(6,5) * t515;
t95 = t243 * t586 + t245 * t527 + t247 * t528 - t360 * t386 + t361 * t388 + t384 * t515;
t856 = t11 - t208 * t513 + t209 * t515 + t586 * t95 + t588 * t94 + (t139 * t621 + t240 * t733) * t617;
t148 = t309 * t525 + t311 * t453 + t313 * t454;
t149 = t310 * t525 + t312 * t453 + t314 * t454;
t187 = t366 * t525 + t367 * t453 + t368 * t454;
t50 = t171 * t525 + t173 * t453 + t175 * t454 + t274 * t311 + t275 * t313 + t309 * t362;
t51 = t170 * t525 + t172 * t453 + t174 * t454 + t274 * t312 + t275 * t314 + t310 * t362;
t81 = t251 * t525 + t252 * t453 + t253 * t454 + t274 * t367 + t275 * t368 + t362 * t366;
t12 = -t148 * t513 + t149 * t515 + t50 * t588 + t51 * t586 + (t187 * t733 + t621 * t81) * t617;
t140 = -t362 * t462 + t363 * t463 + t369 * t588 - t370 * t525 + t371 * t526 - t461 * t513;
t206 = t383 * t588 - t385 * t525 + t387 * t526;
t207 = t384 * t588 - t386 * t525 + t388 * t526;
t239 = t461 * t588 - t462 * t525 + t463 * t526;
t96 = t244 * t588 - t246 * t525 + t248 * t526 - t362 * t385 + t363 * t387 - t383 * t513;
t97 = t243 * t588 - t245 * t525 + t247 * t526 - t362 * t386 + t363 * t388 - t384 * t513;
t855 = t12 - t206 * t513 + t207 * t515 + t586 * t97 + t588 * t96 + (t140 * t621 + t239 * t733) * t617;
t15 = t618 * t80 + (t48 * t622 - t49 * t626 + (t150 * t626 + t151 * t622) * qJD(1)) * t617;
t854 = t15 + t139 * t618 + (t622 * t94 - t626 * t95 + (t208 * t626 + t209 * t622) * qJD(1)) * t617;
t16 = t618 * t81 + (t50 * t622 - t51 * t626 + (t148 * t626 + t149 * t622) * qJD(1)) * t617;
t853 = t16 + t140 * t618 + (t622 * t96 - t626 * t97 + (t206 * t626 + t207 * t622) * qJD(1)) * t617;
t163 = t310 * t561 + t312 * t522 + t314 * t523;
t789 = t163 * t515;
t162 = t309 * t561 + t311 * t522 + t313 * t523;
t790 = t162 * t513;
t100 = t561 * t251 + t522 * t252 + t523 * t253 + t501 * t366 + t378 * t367 + t379 * t368;
t201 = t366 * t561 + t367 * t522 + t368 * t523;
t797 = t100 * t782 + t201 * t705;
t55 = t170 * t561 + t172 * t522 + t174 * t523 + t310 * t501 + t312 * t378 + t314 * t379;
t799 = t55 * t586;
t54 = t171 * t561 + t173 * t522 + t175 * t523 + t309 * t501 + t311 * t378 + t313 * t379;
t800 = t54 * t588;
t22 = t789 - t790 + t797 + t799 + t800;
t159 = t369 * t782 - t561 * t370 + t562 * t371 + t461 * t705 - t501 * t462 + t500 * t463;
t268 = t461 * t782 - t462 * t561 + t463 * t562;
t846 = t268 * t705;
t771 = t159 * t782 + t846;
t218 = t384 * t782 - t386 * t561 + t388 * t562;
t787 = t218 * t515;
t217 = t383 * t782 - t385 * t561 + t387 * t562;
t788 = t217 * t513;
t110 = -t245 * t561 + t247 * t562 - t386 * t501 + t388 * t500 + (t243 * t621 + t384 * t733) * t617;
t791 = t110 * t586;
t109 = -t246 * t561 + t248 * t562 - t385 * t501 + t387 * t500 + (t244 * t621 + t383 * t733) * t617;
t792 = t109 * t588;
t852 = t22 + t771 + t787 - t788 + t791 + t792;
t158 = t159 * t618;
t99 = t100 * t618;
t24 = t99 + (t54 * t622 - t55 * t626 + (t162 * t626 + t163 * t622) * qJD(1)) * t617;
t851 = t24 + t158 + (t109 * t622 - t110 * t626 + (t217 * t626 + t218 * t622) * qJD(1)) * t617;
t73 = t148 * t588 + t149 * t586 + t187 * t782;
t850 = t206 * t588 + t207 * t586 + t239 * t782 + t73;
t74 = t150 * t588 + t151 * t586 + t188 * t782;
t849 = t208 * t588 + t209 * t586 + t240 * t782 + t74;
t76 = t187 * t618 + (t148 * t622 - t149 * t626) * t617;
t848 = t76 + t239 * t618 + (t206 * t622 - t207 * t626) * t617;
t77 = t188 * t618 + (t150 * t622 - t151 * t626) * t617;
t847 = t77 + t240 * t618 + (t208 * t622 - t209 * t626) * t617;
t658 = -t273 * rSges(7,1) - t272 * rSges(7,2);
t176 = t360 * rSges(7,3) - t658;
t807 = t361 * pkin(5);
t770 = t360 * pkin(11) + t176 + t807;
t177 = t275 * rSges(7,1) + t274 * rSges(7,2) + t362 * rSges(7,3);
t769 = t363 * pkin(5) + pkin(11) * t362 + t177;
t254 = rSges(7,1) * t379 + rSges(7,2) * t378 + rSges(7,3) * t501;
t845 = pkin(5) * t500 + pkin(11) * t501 + t254;
t317 = t454 * rSges(7,1) + t453 * rSges(7,2) + t525 * rSges(7,3);
t760 = t526 * pkin(5) + pkin(11) * t525 + t317;
t657 = -t456 * rSges(7,1) - t455 * rSges(7,2);
t318 = -t527 * rSges(7,3) - t657;
t806 = t528 * pkin(5);
t759 = -t527 * pkin(11) + t318 + t806;
t374 = rSges(7,1) * t523 + rSges(7,2) * t522 + rSges(7,3) * t561;
t755 = pkin(5) * t562 + pkin(11) * t561 + t374;
t620 = sin(qJ(4));
t624 = cos(qJ(4));
t538 = t585 * t624 + t620 * t779;
t632 = qJD(4) * t538 + t514 * t620;
t844 = pkin(4) * t632;
t843 = -t830 * t587 - t828 * t588 + t781 * t833;
t842 = -t585 * t831 - t586 * t829 + t779 * t832;
t841 = t830 * t585 + t828 * t586 + t779 * t833;
t840 = t587 * t831 + t588 * t829 + t781 * t832;
t793 = Icges(4,6) * t625;
t556 = Icges(4,4) * t618 + (-Icges(4,2) * t621 - t793) * t617;
t565 = (Icges(4,3) * t621 - t793) * t735;
t566 = (-Icges(4,2) * t625 + t794) * t735;
t777 = t621 * t566;
t839 = ((-t777 + (-qJD(2) * t556 - t565) * t625) * t617 + t857) * t618;
t396 = Icges(4,1) * t707 + Icges(4,4) * t513 - Icges(4,5) * t512;
t397 = -Icges(3,5) * t513 + Icges(3,6) * t512 + Icges(3,3) * t707;
t838 = t397 + t396;
t392 = Icges(4,5) * t707 + Icges(4,6) * t513 - Icges(4,3) * t512;
t399 = -Icges(3,4) * t513 + Icges(3,2) * t512 + Icges(3,6) * t707;
t837 = t399 - t392;
t391 = Icges(4,5) * t708 - Icges(4,6) * t515 + Icges(4,3) * t514;
t400 = Icges(3,4) * t515 - Icges(3,2) * t514 + Icges(3,6) * t708;
t836 = t400 - t391;
t394 = Icges(4,4) * t707 + Icges(4,2) * t513 - Icges(4,6) * t512;
t401 = -Icges(3,1) * t513 + Icges(3,4) * t512 + Icges(3,5) * t707;
t835 = t401 - t394;
t393 = Icges(4,4) * t708 - Icges(4,2) * t515 + Icges(4,6) * t514;
t402 = Icges(3,1) * t515 - Icges(3,4) * t514 + Icges(3,5) * t708;
t834 = t402 - t393;
t395 = Icges(4,1) * t708 - Icges(4,4) * t515 + Icges(4,5) * t514;
t398 = Icges(3,5) * t515 - Icges(3,6) * t514 + Icges(3,3) * t708;
t827 = (-t398 - t395) * t626;
t616 = t617 ^ 2;
t826 = m(7) / 0.2e1;
t825 = t360 / 0.2e1;
t824 = t362 / 0.2e1;
t823 = t501 / 0.2e1;
t822 = -t513 / 0.2e1;
t821 = t515 / 0.2e1;
t820 = t525 / 0.2e1;
t819 = -t527 / 0.2e1;
t818 = t561 / 0.2e1;
t817 = t586 / 0.2e1;
t816 = t588 / 0.2e1;
t815 = t617 / 0.2e1;
t814 = t618 / 0.2e1;
t813 = t622 / 0.2e1;
t812 = -t626 / 0.2e1;
t811 = t626 / 0.2e1;
t810 = rSges(4,2) - pkin(2);
t809 = -rSges(5,3) - pkin(2);
t808 = rSges(7,3) + pkin(11);
t615 = t626 * pkin(1);
t627 = -pkin(10) - pkin(9);
t805 = -pkin(2) + t627;
t613 = pkin(4) * t624 + pkin(3);
t804 = -pkin(3) + t613;
t803 = -pkin(9) - t627;
t802 = rSges(4,3) * t514;
t801 = pkin(4) * qJD(4);
t798 = t100 * t561 + t201 * t501;
t786 = t512 * t620;
t784 = t585 * t620;
t783 = t587 * t620;
t778 = t620 * t625;
t583 = -t618 * t620 - t624 * t780;
t533 = qJD(4) * t583 + t620 * t704;
t648 = t617 * t778 - t618 * t624;
t534 = qJD(4) * t648 + t624 * t704;
t440 = Icges(5,5) * t533 + Icges(5,6) * t534 + Icges(5,3) * t705;
t441 = Icges(5,4) * t533 + Icges(5,2) * t534 + Icges(5,6) * t705;
t442 = Icges(5,1) * t533 + Icges(5,4) * t534 + Icges(5,5) * t705;
t471 = -Icges(5,5) * t648 + Icges(5,6) * t583 + Icges(5,3) * t782;
t472 = -Icges(5,4) * t648 + Icges(5,2) * t583 + Icges(5,6) * t782;
t473 = -Icges(5,1) * t648 + Icges(5,4) * t583 + Icges(5,5) * t782;
t186 = t440 * t782 + t583 * t441 - t442 * t648 + t471 * t705 + t534 * t472 + t533 * t473;
t295 = t471 * t782 + t472 * t583 - t473 * t648;
t768 = t186 * t782 + t295 * t705;
t660 = -t361 * rSges(6,1) + t360 * rSges(6,2);
t249 = t515 * rSges(6,3) - t660;
t659 = -t528 * rSges(6,1) - t527 * rSges(6,2);
t390 = t586 * rSges(6,3) - t659;
t767 = t588 * t249 - t513 * t390;
t250 = t363 * rSges(6,1) - t362 * rSges(6,2) - t513 * rSges(6,3);
t389 = t526 * rSges(6,1) - t525 * rSges(6,2) + t588 * rSges(6,3);
t766 = t250 * t782 + t389 * t705;
t511 = t515 * pkin(9);
t293 = -t515 * t627 + t708 * t804 - t511 + t844;
t667 = -pkin(4) * t784 + t613 * t779;
t738 = pkin(3) * t779 - t586 * pkin(9);
t434 = -t586 * t627 - t667 + t738;
t764 = t588 * t293 - t513 * t434;
t726 = t624 * t801;
t671 = -pkin(4) * t786 + t513 * t627 + t587 * t726 + t613 * t707;
t727 = t620 * t801;
t680 = t617 * t727;
t691 = -pkin(3) * t707 + pkin(9) * t513;
t294 = -t622 * t680 + t671 + t691;
t543 = pkin(3) * t781 + pkin(9) * t588;
t709 = pkin(4) * t783 - t588 * t627 + t613 * t781;
t433 = -t543 + t709;
t763 = t294 * t782 + t433 * t705;
t762 = t759 * t588;
t761 = t760 * t782;
t758 = t755 * t586;
t375 = rSges(6,1) * t500 - rSges(6,2) * t501 + rSges(6,3) * t705;
t465 = rSges(6,1) * t562 - rSges(6,2) * t561 + rSges(6,3) * t782;
t757 = t586 * t375 + t515 * t465;
t746 = -t514 * qJ(3) - t585 * qJD(3);
t358 = t515 * pkin(2) - t746;
t469 = pkin(3) * t708 + t511;
t756 = -t358 - t469;
t467 = -t618 * t727 + (-t625 * t726 + (pkin(4) * t620 * t621 + t625 * t803) * qJD(2)) * t617;
t754 = -t375 - t467;
t753 = -t389 - t433;
t752 = -t390 - t434;
t497 = t804 * t618 + (-pkin(4) * t778 + t621 * t803) * t617;
t751 = t586 * t467 + t515 * t497;
t443 = rSges(5,1) * t533 + rSges(5,2) * t534 + rSges(5,3) * t705;
t542 = (-qJD(3) * t625 + (pkin(2) * t625 + qJ(3) * t621) * qJD(2)) * t617;
t750 = -t443 - t542;
t749 = t465 + t497;
t574 = t585 * qJ(3);
t516 = t586 * pkin(2) + t574;
t517 = t588 * pkin(2) + qJ(3) * t587;
t748 = t516 * t781 + t517 * t779;
t502 = t618 * t517;
t747 = t618 * t543 + t502;
t745 = -t516 + t738;
t744 = -t517 - t543;
t742 = -t542 - (-rSges(4,2) * t625 + rSges(4,3) * t621) * t735;
t559 = rSges(4,1) * t618 + (-rSges(4,2) * t621 - rSges(4,3) * t625) * t617;
t589 = (pkin(2) * t621 - qJ(3) * t625) * t617;
t741 = -t559 - t589;
t593 = t618 * pkin(3) + pkin(9) * t782;
t740 = -t589 - t593;
t739 = pkin(8) * t781 + t615;
t729 = -rSges(6,3) + t805;
t728 = pkin(1) * t737;
t725 = t54 / 0.2e1 + t81 / 0.2e1;
t724 = -t55 / 0.2e1 - t80 / 0.2e1;
t721 = -t467 - t845;
t720 = -t293 + t756;
t719 = -t433 - t760;
t718 = -t434 - t759;
t357 = -t513 * pkin(2) - qJ(3) * t512 + qJD(3) * t587;
t717 = t357 * t779 + t358 * t781 + t516 * t707;
t716 = t497 + t755;
t715 = -t542 + t754;
t537 = t624 * t781 + t783;
t418 = -qJD(4) * t537 - t512 * t624 - t620 * t707;
t536 = t587 * t624 - t620 * t781;
t419 = qJD(4) * t536 + t624 * t707 - t786;
t267 = t419 * rSges(5,1) + t418 * rSges(5,2) - t513 * rSges(5,3);
t714 = t618 * t433 + t747;
t713 = -t433 + t744;
t712 = -t434 + t745;
t474 = -rSges(5,1) * t648 + rSges(5,2) * t583 + rSges(5,3) * t782;
t711 = -t474 + t740;
t710 = -t497 + t740;
t405 = -t513 * rSges(3,1) + t512 * rSges(3,2) + rSges(3,3) * t707;
t428 = t537 * rSges(5,1) + t536 * rSges(5,2) + t588 * rSges(5,3);
t492 = t588 * rSges(3,1) - t587 * rSges(3,2) + rSges(3,3) * t781;
t404 = rSges(4,1) * t707 + t513 * rSges(4,2) - t512 * rSges(4,3);
t489 = rSges(4,1) * t781 - t588 * rSges(4,2) + t587 * rSges(4,3);
t706 = t616 * t733;
t703 = t782 / 0.2e1;
t261 = Icges(5,5) * t419 + Icges(5,6) * t418 - Icges(5,3) * t513;
t263 = Icges(5,4) * t419 + Icges(5,2) * t418 - Icges(5,6) * t513;
t265 = Icges(5,1) * t419 + Icges(5,4) * t418 - Icges(5,5) * t513;
t422 = Icges(5,5) * t537 + Icges(5,6) * t536 + Icges(5,3) * t588;
t424 = Icges(5,4) * t537 + Icges(5,2) * t536 + Icges(5,6) * t588;
t426 = Icges(5,1) * t537 + Icges(5,4) * t536 + Icges(5,5) * t588;
t120 = t263 * t583 - t265 * t648 + t424 * t534 + t426 * t533 + (t261 * t621 + t422 * t733) * t617;
t145 = t418 * t472 + t419 * t473 + t440 * t588 + t441 * t536 + t442 * t537 - t471 * t513;
t700 = t120 / 0.2e1 + t145 / 0.2e1;
t649 = t624 * t779 - t784;
t416 = qJD(4) * t649 + t514 * t624 - t620 * t708;
t417 = t624 * t708 + t632;
t260 = Icges(5,5) * t417 + Icges(5,6) * t416 + Icges(5,3) * t515;
t262 = Icges(5,4) * t417 + Icges(5,2) * t416 + Icges(5,6) * t515;
t264 = Icges(5,1) * t417 + Icges(5,4) * t416 + Icges(5,5) * t515;
t423 = -Icges(5,5) * t649 + Icges(5,6) * t538 + Icges(5,3) * t586;
t425 = -Icges(5,4) * t649 + Icges(5,2) * t538 + Icges(5,6) * t586;
t427 = -Icges(5,1) * t649 + Icges(5,4) * t538 + Icges(5,5) * t586;
t121 = t262 * t583 - t264 * t648 + t425 * t534 + t427 * t533 + (t260 * t621 + t423 * t733) * t617;
t144 = t416 * t472 + t417 * t473 + t440 * t586 + t441 * t538 - t442 * t649 + t471 * t515;
t699 = t121 / 0.2e1 + t144 / 0.2e1;
t698 = t162 / 0.2e1 + t187 / 0.2e1;
t697 = t163 / 0.2e1 + t188 / 0.2e1;
t226 = t423 * t782 + t425 * t583 - t427 * t648;
t257 = t471 * t586 + t472 * t538 - t473 * t649;
t696 = t226 / 0.2e1 + t257 / 0.2e1;
t225 = t422 * t782 + t424 * t583 - t426 * t648;
t256 = t471 * t588 + t472 * t536 + t473 * t537;
t695 = -t256 / 0.2e1 - t225 / 0.2e1;
t692 = -t622 * pkin(1) + pkin(8) * t779;
t690 = 2 * m(3);
t689 = 2 * m(4);
t687 = 2 * m(5);
t685 = 2 * m(6);
t683 = 0.2e1 * m(7);
t682 = t741 * t626;
t679 = pkin(9) * t706;
t677 = -t513 * t759 + t588 * t770;
t676 = t705 * t760 + t769 * t782;
t675 = t515 * t755 + t586 * t845;
t674 = -t542 + t721;
t673 = -t465 + t710;
t672 = t543 * t779 - t738 * t781 + t748;
t668 = t705 / 0.2e1;
t666 = -t574 + t692;
t665 = t711 * t626;
t663 = -rSges(3,1) * t515 + rSges(3,2) * t514;
t662 = -rSges(5,1) * t417 - rSges(5,2) * t416;
t661 = rSges(5,1) * t649 - rSges(5,2) * t538;
t656 = t710 - t755;
t655 = t673 * t626;
t654 = t517 + t739;
t652 = rSges(4,1) * t779 - rSges(4,3) * t585;
t651 = t469 * t781 - t691 * t779 - t707 * t738 + t717;
t650 = t433 * t779 + t434 * t781 + t672;
t646 = t656 * t626;
t603 = pkin(8) * t707;
t645 = t357 + t603;
t352 = t618 * t357;
t644 = -t618 * t691 - t622 * t679 + t352;
t547 = t589 * t708;
t643 = t593 * t708 - t626 * t679 + t547;
t491 = rSges(3,1) * t586 - rSges(3,2) * t585 - rSges(3,3) * t779;
t641 = t618 * t294 + t644;
t640 = t497 * t708 + t643;
t639 = t666 + t667;
t638 = t654 + t709;
t636 = t293 * t781 + t294 * t779 + t434 * t707 + t651;
t18 = t162 * t362 + t163 * t360 + t54 * t525 - t55 * t527 + t798;
t3 = t150 * t362 + t151 * t360 + t188 * t501 + t48 * t525 - t49 * t527 + t561 * t80;
t4 = t148 * t362 + t149 * t360 + t187 * t501 + t50 * t525 - t51 * t527 + t561 * t81;
t63 = t148 * t525 - t149 * t527 + t187 * t561;
t64 = t150 * t525 - t151 * t527 + t188 * t561;
t84 = t162 * t525 - t163 * t527 + t201 * t561;
t87 = t162 * t588 + t163 * t586 + t201 * t782;
t635 = t11 * t819 + t12 * t820 + t18 * t703 + t22 * t818 + t3 * t817 + t4 * t816 + t63 * t822 + t64 * t821 + t84 * t668 + t73 * t824 + t74 * t825 + t87 * t823;
t634 = t87 * t705 + (t846 + t852) * t782 + (t217 * t705 + t855) * t588 + (t218 * t705 + t856) * t586 + t849 * t515 - t850 * t513;
t633 = t645 - t728;
t631 = t792 / 0.2e1 + t791 / 0.2e1 - t790 / 0.2e1 + t789 / 0.2e1 - t788 / 0.2e1 + t787 / 0.2e1 + t800 / 0.2e1 + t799 / 0.2e1 + t771 + t797 + (t187 + t239) * t822 + (t188 + t240) * t821 + (t139 + t80) * t817 + (t140 + t81) * t816;
t91 = t201 * t618 + (t162 * t622 - t163 * t626) * t617;
t630 = t848 * t822 + t847 * t821 + t854 * t817 + t853 * t816 + t852 * t814 + t851 * t703 + t855 * t781 / 0.2e1 - t856 * t779 / 0.2e1 + (t268 * t618 + (t217 * t622 - t218 * t626) * t617 + t91) * t668 + (t622 * t849 + t626 * t850) * qJD(1) * t815;
t629 = (-pkin(1) * qJD(1) - t680) * t622 + t645 + t671;
t628 = -t844 + (-t615 + (-pkin(8) - t613) * t781) * qJD(1) + t746;
t572 = (rSges(3,1) * t625 - rSges(3,2) * t621) * t735;
t558 = rSges(3,3) * t618 + (rSges(3,1) * t621 + rSges(3,2) * t625) * t617;
t557 = Icges(4,1) * t618 + (-Icges(4,4) * t621 - Icges(4,5) * t625) * t617;
t552 = Icges(3,3) * t618 + (Icges(3,5) * t621 + Icges(3,6) * t625) * t617;
t490 = -rSges(4,2) * t586 - t652;
t458 = t492 + t739;
t457 = -t491 + t692;
t452 = t586 * t497;
t446 = t586 * t465;
t439 = -t491 * t618 - t558 * t779;
t438 = t492 * t618 - t558 * t781;
t429 = rSges(5,3) * t586 - t661;
t414 = t433 * t782;
t406 = rSges(3,3) * t708 - t663;
t403 = rSges(4,1) * t708 - rSges(4,2) * t515 + t802;
t377 = t654 + t489;
t376 = t586 * t810 + t652 + t666;
t365 = t588 * t434;
t364 = t389 * t782;
t349 = (-t615 + (-rSges(3,3) - pkin(8)) * t781) * qJD(1) + t663;
t348 = t405 + t603 - t728;
t345 = t552 * t781 - t553 * t587 + t554 * t588;
t344 = -t552 * t779 - t553 * t585 + t554 * t586;
t343 = t555 * t585 - t556 * t586 - t557 * t779;
t342 = t555 * t587 - t556 * t588 + t557 * t781;
t341 = t588 * t390;
t328 = (-t490 - t516) * t618 + t617 * t682;
t327 = t489 * t618 + t741 * t781 + t502;
t326 = t405 * t618 + (-t558 * t736 - t572 * t622) * t617;
t325 = -t406 * t618 + (t558 * t737 - t572 * t626) * t617;
t324 = t428 * t782 - t474 * t588;
t323 = -t429 * t782 + t474 * t586;
t322 = t481 * t618 + (-t477 * t625 - t479 * t621) * t617;
t321 = t480 * t618 + (-t476 * t625 - t478 * t621) * t617;
t320 = t483 * t618 + (t485 * t625 + t487 * t621) * t617;
t319 = t482 * t618 + (t484 * t625 + t486 * t621) * t617;
t316 = t543 + t654 + t428;
t315 = t586 * t809 + t661 + t666 + t738;
t306 = -t465 * t588 + t364;
t305 = -t390 * t782 + t446;
t285 = (t489 * t626 + t490 * t622) * t617 + t748;
t284 = t638 + t389;
t283 = t586 * t729 + t639 + t659;
t278 = -t428 * t586 + t429 * t588;
t277 = -t802 + t810 * t515 + (-t615 + (-rSges(4,1) - pkin(8)) * t781) * qJD(1) + t746;
t276 = t633 + t404;
t266 = rSges(5,3) * t515 - t662;
t259 = -t389 * t586 + t341;
t238 = -t514 * t553 + t515 * t554 - t569 * t585 + t570 * t586 + (t552 * t737 - t567 * t626) * t617;
t237 = t512 * t553 - t513 * t554 - t569 * t587 + t570 * t588 + (t552 * t736 + t567 * t622) * t617;
t236 = -t512 * t555 + t513 * t556 + t565 * t587 - t566 * t588 + (t557 * t736 + t568 * t622) * t617;
t235 = t514 * t555 - t515 * t556 + t565 * t585 - t566 * t586 + (t557 * t737 - t568 * t626) * t617;
t234 = (-t429 + t745) * t618 + t617 * t665;
t233 = t428 * t618 + t711 * t781 + t747;
t224 = -t318 * t561 - t374 * t527;
t223 = t317 * t561 - t374 * t525;
t222 = t404 * t618 + t352 + (qJD(1) * t682 + t622 * t742) * t617;
t221 = t547 + (-t358 - t403) * t618 + (t559 * t737 + t626 * t742) * t617;
t220 = -t588 * t749 + t364 + t414;
t219 = t752 * t782 + t446 + t452;
t216 = t423 * t586 + t425 * t538 - t427 * t649;
t215 = t422 * t586 + t424 * t538 - t426 * t649;
t214 = t423 * t588 + t425 * t536 + t427 * t537;
t213 = t422 * t588 + t424 * t536 + t426 * t537;
t212 = t638 + t760;
t211 = t527 * t808 + t586 * t805 + t639 + t657 - t806;
t210 = (t428 * t626 + t429 * t622) * t617 + t672;
t203 = -t511 + t809 * t515 + (-t615 + (-pkin(3) - pkin(8)) * t781) * qJD(1) + t662 + t746;
t202 = t633 - t691 + t267;
t200 = t317 * t527 + t318 * t525;
t198 = t396 * t618 + (-t392 * t625 - t394 * t621 + (t476 * t621 - t478 * t625) * qJD(2)) * t617;
t197 = t395 * t618 + (-t391 * t625 - t393 * t621 + (t477 * t621 - t479 * t625) * qJD(2)) * t617;
t196 = t397 * t618 + (t399 * t625 + t401 * t621 + (-t485 * t621 + t487 * t625) * qJD(2)) * t617;
t195 = t398 * t618 + (t400 * t625 + t402 * t621 + (-t484 * t621 + t486 * t625) * qJD(2)) * t617;
t194 = (-t390 + t712) * t618 + t617 * t655;
t193 = t389 * t618 + t673 * t781 + t714;
t192 = t586 * t753 + t341 + t365;
t191 = -t588 * t755 + t761;
t190 = -t759 * t782 + t758;
t185 = t186 * t618;
t181 = t515 * t729 + t628 + t660;
t180 = t629 + t250;
t179 = t443 * t586 + t474 * t515 + (-t266 * t621 - t429 * t733) * t617;
t178 = -t443 * t588 + t474 * t513 + (t267 * t621 + t428 * t733) * t617;
t168 = (t389 * t626 + t390 * t622) * t617 + t650;
t167 = -t586 * t760 + t762;
t165 = (-t249 * t621 - t390 * t733) * t617 + t757;
t164 = -t375 * t588 + t465 * t513 + t766;
t161 = -t588 * t716 + t414 + t761;
t160 = t718 * t782 + t452 + t758;
t156 = (t403 * t622 + t404 * t626 + (t490 * t626 + (-t489 - t517) * t622) * qJD(1)) * t617 + t717;
t155 = t267 * t618 + (qJD(1) * t665 + t622 * t750) * t617 + t644;
t154 = (-t266 + t756) * t618 + (t474 * t737 + t626 * t750) * t617 + t643;
t147 = (t712 - t759) * t618 + t617 * t646;
t146 = t618 * t760 + t656 * t781 + t714;
t143 = t266 * t588 - t267 * t586 - t428 * t515 - t429 * t513;
t142 = t586 * t719 + t365 + t762;
t141 = -t250 * t586 - t389 * t515 + t767;
t136 = (t622 * t759 + t626 * t760) * t617 + t650;
t134 = t257 * t618 + (t215 * t622 - t216 * t626) * t617;
t133 = t256 * t618 + (t213 * t622 - t214 * t626) * t617;
t131 = t215 * t588 + t216 * t586 + t257 * t782;
t130 = t213 * t588 + t214 * t586 + t256 * t782;
t129 = t250 * t618 + (qJD(1) * t655 + t622 * t715) * t617 + t641;
t128 = (-t249 + t720) * t618 + (t465 * t737 + t626 * t715) * t617 + t640;
t127 = -t360 * t808 + t515 * t805 + t628 + t658 - t807;
t126 = t629 + t769;
t125 = ((-t249 - t293) * t621 + t752 * t733) * t617 + t751 + t757;
t124 = t513 * t749 + t588 * t754 + t763 + t766;
t111 = (t266 * t622 + t267 * t626 + (t429 * t626 + (-t428 + t744) * t622) * qJD(1)) * t617 + t651;
t106 = t260 * t588 + t262 * t536 + t264 * t537 + t418 * t425 + t419 * t427 - t423 * t513;
t105 = t261 * t588 + t263 * t536 + t265 * t537 + t418 * t424 + t419 * t426 - t422 * t513;
t104 = t260 * t586 + t262 * t538 - t264 * t649 + t416 * t425 + t417 * t427 + t423 * t515;
t103 = t261 * t586 + t263 * t538 - t265 * t649 + t416 * t424 + t417 * t426 + t422 * t515;
t102 = t177 * t561 - t254 * t525 + t317 * t501 - t362 * t374;
t101 = -t176 * t561 - t254 * t527 - t318 * t501 + t360 * t374;
t92 = (-t250 - t294) * t586 + t753 * t515 + t764 + t767;
t89 = (-t621 * t770 - t733 * t759) * t617 + t675;
t88 = t513 * t755 - t588 * t845 + t676;
t83 = t176 * t525 + t177 * t527 - t317 * t360 + t318 * t362;
t75 = (t249 * t622 + t250 * t626 + (t390 * t626 + (-t389 + t713) * t622) * qJD(1)) * t617 + t636;
t68 = t769 * t618 + (qJD(1) * t646 + t622 * t674) * t617 + t641;
t67 = (t720 - t770) * t618 + (t626 * t674 + t737 * t755) * t617 + t640;
t58 = ((-t293 - t770) * t621 + t718 * t733) * t617 + t675 + t751;
t57 = t513 * t716 + t588 * t721 + t676 + t763;
t56 = -t515 * t760 - t586 * t769 + t677;
t47 = (-t294 - t769) * t586 + t719 * t515 + t677 + t764;
t46 = (t769 * t626 + t770 * t622 + (t759 * t626 + (t713 - t760) * t622) * qJD(1)) * t617 + t636;
t45 = t185 + (t120 * t622 - t121 * t626 + (t225 * t626 + t226 * t622) * qJD(1)) * t617;
t44 = t120 * t588 + t121 * t586 - t225 * t513 + t226 * t515 + t768;
t41 = t145 * t618 + (t105 * t622 - t106 * t626 + (t213 * t626 + t214 * t622) * qJD(1)) * t617;
t40 = t144 * t618 + (t103 * t622 - t104 * t626 + (t215 * t626 + t216 * t622) * qJD(1)) * t617;
t36 = t105 * t588 + t106 * t586 - t213 * t513 + t214 * t515 + (t145 * t621 + t256 * t733) * t617;
t35 = t103 * t588 + t104 * t586 - t215 * t513 + t216 * t515 + (t144 * t621 + t257 * t733) * t617;
t1 = [t186 + (t276 * t377 + t277 * t376) * t689 + (t348 * t458 + t349 * t457) * t690 - t565 * t780 - t556 * t705 - t617 * t777 + (t126 * t212 + t127 * t211) * t683 + (t180 * t284 + t181 * t283) * t685 + (t202 * t316 + t203 * t315) * t687 + t159 + t100 + t857; (t154 * t315 + t155 * t316 + t202 * t233 + t203 * t234) * m(5) + (t128 * t283 + t129 * t284 + t180 * t193 + t181 * t194) * m(6) + (t126 * t146 + t127 * t147 + t211 * t67 + t212 * t68) * m(7) + t99 + t158 + t185 + ((-t195 / 0.2e1 - t197 / 0.2e1 - t139 / 0.2e1 - t110 / 0.2e1 - t238 / 0.2e1 - t235 / 0.2e1 - t699 + t724) * t626 + (t196 / 0.2e1 + t198 / 0.2e1 + t109 / 0.2e1 + t140 / 0.2e1 + t237 / 0.2e1 + t236 / 0.2e1 + t700 + t725) * t622 + ((t320 / 0.2e1 + t321 / 0.2e1 + t239 / 0.2e1 + t217 / 0.2e1 + t342 / 0.2e1 + t345 / 0.2e1 - t695 + t698) * t626 + (t319 / 0.2e1 + t322 / 0.2e1 + t240 / 0.2e1 + t218 / 0.2e1 + t343 / 0.2e1 + t344 / 0.2e1 + t696 + t697) * t622) * qJD(1)) * t617 + m(3) * (t325 * t457 + t326 * t458 + t348 * t438 + t349 * t439) + m(4) * (t221 * t376 + t222 * t377 + t276 * t327 + t277 * t328) + t839; (t111 * t210 + t154 * t234 + t155 * t233) * t687 + (t156 * t285 + t221 * t328 + t222 * t327) * t689 + (t128 * t194 + t129 * t193 + t168 * t75) * t685 + (t439 * t325 + t438 * t326 + (t491 * t622 + t492 * t626) * (t405 * t626 + t406 * t622 + (t491 * t626 - t492 * t622) * qJD(1)) * t616) * t690 + (t136 * t46 + t146 * t68 + t147 * t67) * t683 + (t41 + t853) * t781 + (-t40 - t854) * t779 + (t134 + t847) * t708 + (t133 + t848) * t707 + ((-t622 * t841 + t626 * t842) * t708 + (t622 * t843 - t626 * t840) * t707 + (t840 * t737 + t843 * t736 + (t512 * t831 + t513 * t829 + t587 * t836 - t588 * t834 - t707 * t832) * t626 + ((t622 * t838 + t736 * t833 + t827) * t617 + t835 * t588 - t837 * t587 + t828 * t513 + t830 * t512) * t622) * t781 + (t842 * t737 + t841 * t736 + ((t737 * t832 + t827) * t617 + t834 * t586 - t836 * t585 + t829 * t515 + t831 * t514) * t626 + ((t626 * t838 - t737 * t833) * t617 - t835 * t586 + t837 * t585 + t828 * t515 + t830 * t514) * t622) * t779) * t617 + (t45 + (t237 + t236) * t781 + (-t238 - t235) * t779 + (t344 + t343) * t708 + (t345 + t342) * t707 + ((-t195 - t197) * t626 + (t196 + t198) * t622 + ((t320 + t321) * t626 + (t319 + t322) * t622) * qJD(1)) * t617 + t839 + t851) * t618; (m(4) * t277 + m(5) * t203 + m(6) * t181 + m(7) * t127) * t587 + (m(4) * t276 + m(5) * t202 + m(6) * t180 + m(7) * t126) * t585 + (m(4) * t377 + m(5) * t316 + m(6) * t284 + m(7) * t212) * t514 + (-m(4) * t376 - m(5) * t315 - m(6) * t283 - m(7) * t211) * t512; (-m(4) * t156 - m(5) * t111 - m(6) * t75 - m(7) * t46) * t780 + (m(4) * t221 + m(5) * t154 + m(6) * t128 + m(7) * t67) * t587 + (m(4) * t222 + m(5) * t155 + m(6) * t129 + m(7) * t68) * t585 + (m(4) * t327 + m(5) * t233 + m(6) * t193 + m(7) * t146) * t514 + (-m(4) * t328 - m(5) * t234 - m(6) * t194 - m(7) * t147) * t512 + (m(4) * t285 + m(5) * t210 + m(6) * t168 + m(7) * t136) * t704; 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + t826) * (-t512 * t587 + t514 * t585 - t621 * t706); t631 + t699 * t586 + t696 * t515 + (t124 * t284 + t125 * t283 + t180 * t220 + t181 * t219) * m(6) + (t126 * t161 + t127 * t160 + t211 * t58 + t212 * t57) * m(7) + (t178 * t316 + t179 * t315 + t202 * t324 + t203 * t323) * m(5) + t695 * t513 + t700 * t588 + t768; (t124 * t193 + t125 * t194 + t128 * t219 + t129 * t220 + t168 * t92 + t192 * t75) * m(6) + (t136 * t47 + t142 * t46 + t146 * t57 + t147 * t58 + t160 * t67 + t161 * t68) * m(7) + t134 * t821 + t133 * t822 + t630 + (t111 * t278 + t143 * t210 + t154 * t323 + t155 * t324 + t178 * t233 + t179 * t234) * m(5) + (t35 * t812 + t36 * t813 + t621 * t45 / 0.2e1 + (t295 * t814 + (t225 * t622 - t226 * t626) * t815) * t733 + (t130 * t811 + t131 * t813) * qJD(1)) * t617 + t40 * t817 + t41 * t816 + t44 * t814; (-m(5) * t143 - m(6) * t92 - m(7) * t47) * t780 + (m(5) * t179 + m(6) * t125 + m(7) * t58) * t587 + (m(5) * t178 + m(6) * t124 + m(7) * t57) * t585 + (m(5) * t324 + m(6) * t220 + m(7) * t161) * t514 + (-m(5) * t323 - m(6) * t219 - m(7) * t160) * t512 + (m(5) * t278 + m(6) * t192 + m(7) * t142) * t704; t515 * t131 - t513 * t130 + (t124 * t220 + t125 * t219 + t192 * t92) * t685 + (t142 * t47 + t160 * t58 + t161 * t57) * t683 + (t143 * t278 + t178 * t324 + t179 * t323) * t687 + (t225 * t588 + t226 * t586 + t295 * t782) * t705 + t44 * t782 + t586 * t35 + t634 + t588 * t36; t631 + (t126 * t191 + t127 * t190 + t211 * t89 + t212 * t88) * m(7) + (t164 * t284 + t165 * t283 + t180 * t306 + t181 * t305) * m(6); (t128 * t305 + t129 * t306 + t141 * t168 + t164 * t193 + t165 * t194 + t259 * t75) * m(6) + (t136 * t56 + t146 * t88 + t147 * t89 + t167 * t46 + t190 * t67 + t191 * t68) * m(7) + t630; (-m(6) * t141 - m(7) * t56) * t780 + (m(6) * t165 + m(7) * t89) * t587 + (m(6) * t164 + m(7) * t88) * t585 + (m(6) * t306 + m(7) * t191) * t514 + (-m(6) * t305 - m(7) * t190) * t512 + (m(6) * t259 + m(7) * t167) * t704; (t124 * t306 + t125 * t305 + t141 * t192 + t164 * t220 + t165 * t219 + t259 * t92) * m(6) + (t142 * t56 + t160 * t89 + t161 * t88 + t167 * t47 + t190 * t58 + t191 * t57) * m(7) + t634; (t141 * t259 + t164 * t306 + t165 * t305) * t685 + (t167 * t56 + t190 * t89 + t191 * t88) * t683 + t634; (t101 * t211 + t102 * t212 + t126 * t223 + t127 * t224) * m(7) + t724 * t527 + t725 * t525 + t698 * t362 + t697 * t360 + t798; t18 * t814 + t91 * t823 + t24 * t818 + t76 * t824 + t16 * t820 + t77 * t825 + t15 * t819 + (t101 * t147 + t102 * t146 + t136 * t83 + t200 * t46 + t223 * t68 + t224 * t67) * m(7) + (t3 * t812 + t4 * t813 + (t63 * t811 + t64 * t813) * qJD(1)) * t617; 0.2e1 * (t101 * t587 + t102 * t585 + t223 * t514 - t224 * t512 + (t200 * t734 - t625 * t83) * t617) * t826; t635 + (t101 * t160 + t102 * t161 + t142 * t83 + t200 * t47 + t223 * t57 + t224 * t58) * m(7); t635 + (t101 * t190 + t102 * t191 + t167 * t83 + t200 * t56 + t223 * t88 + t224 * t89) * m(7); t362 * t63 + t525 * t4 + t360 * t64 - t527 * t3 + t501 * t84 + t561 * t18 + (t101 * t224 + t102 * t223 + t200 * t83) * t683;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;