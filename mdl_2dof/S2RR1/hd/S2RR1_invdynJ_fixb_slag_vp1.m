% Calculate vector of inverse dynamics joint torques for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_invdynJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:06
% EndTime: 2020-01-03 11:19:09
% DurationCPUTime: 1.80s
% Computational Cost: add. (1143->213), mult. (3010->323), div. (0->0), fcn. (2394->4), ass. (0->128)
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t84 = cos(qJ(2));
t137 = Icges(3,4) * t84;
t82 = sin(qJ(2));
t98 = -Icges(3,2) * t82 + t137;
t42 = Icges(3,6) * t85 + t98 * t83;
t133 = Icges(3,2) * t84;
t138 = Icges(3,4) * t82;
t68 = -t133 - t138;
t70 = -Icges(3,1) * t82 - t137;
t100 = t82 * t68 - t84 * t70;
t136 = Icges(3,5) * t84;
t97 = -Icges(3,6) * t82 + t136;
t161 = t100 * qJD(1) + t97 * qJD(2);
t132 = Icges(3,6) * t83;
t43 = t98 * t85 - t132;
t150 = t43 * t82;
t139 = Icges(3,1) * t84;
t99 = -t138 + t139;
t45 = -Icges(3,5) * t83 + t99 * t85;
t103 = -t45 * t84 + t150;
t130 = Icges(3,3) * t85;
t66 = -Icges(3,5) * t82 - Icges(3,6) * t84;
t93 = qJD(2) * t66;
t160 = -t85 * t93 + (t97 * t83 + t103 + t130) * qJD(1);
t135 = Icges(3,5) * t85;
t79 = t83 * t138;
t44 = t83 * t139 + t135 - t79;
t104 = t42 * t82 - t44 * t84;
t41 = -Icges(3,3) * t83 + t97 * t85;
t129 = qJD(1) * t41;
t159 = t104 * qJD(1) - t83 * t93 - t129;
t109 = rSges(3,1) * t82 + rSges(3,2) * t84;
t127 = qJD(2) * t85;
t147 = t83 * t84;
t122 = rSges(3,1) * t147;
t149 = t82 * t83;
t80 = rSges(3,2) * t149;
t115 = t80 - t122;
t46 = rSges(3,3) * t85 - t115;
t24 = -(pkin(1) * t85 + t46) * qJD(1) - t109 * t127;
t158 = -t42 * t85 + t43 * t83;
t123 = qJD(1) * qJD(2);
t64 = -qJDD(2) * t83 - t85 * t123;
t157 = t64 / 0.2e1;
t65 = qJDD(2) * t85 - t83 * t123;
t156 = t65 / 0.2e1;
t155 = rSges(3,3) + pkin(1);
t154 = rSges(3,1) * t84;
t153 = rSges(3,3) * t83;
t124 = t83 * qJD(2);
t128 = qJD(1) * t83;
t146 = t84 * t85;
t148 = t82 * t85;
t114 = rSges(3,1) * t146 - rSges(3,2) * t148;
t47 = t114 - t153;
t23 = pkin(1) * t128 - qJD(1) * t47 + t109 * t124;
t55 = t109 * t85;
t152 = t23 * t55;
t151 = t24 * t83;
t48 = t66 * t83;
t49 = t66 * t85;
t38 = t85 * t41;
t144 = t45 * t147 + t38;
t141 = t68 + t99;
t140 = t98 - t70;
t21 = -t100 * t83 + t49;
t126 = t21 * qJD(1);
t125 = t97 * qJD(1);
t120 = -t127 / 0.2e1;
t119 = t127 / 0.2e1;
t118 = -t124 / 0.2e1;
t117 = t124 / 0.2e1;
t40 = -t82 * t132 + t83 * t136 + t130;
t116 = -t40 - t150;
t113 = t44 * t146 - t42 * t148;
t20 = -t43 * t84 - t45 * t82;
t94 = qJD(2) * t68;
t27 = -t42 * qJD(1) + t85 * t94;
t95 = qJD(2) * t70;
t29 = t85 * t95 + (-t99 * t83 - t135) * qJD(1);
t89 = t20 * qJD(2) - t27 * t82 + t29 * t84 - t129;
t19 = -t42 * t84 - t44 * t82;
t28 = t43 * qJD(1) + t83 * t94;
t30 = t45 * qJD(1) + t83 * t95;
t90 = -qJD(1) * t40 + t19 * qJD(2) - t28 * t82 + t30 * t84;
t112 = (t159 * t83 + t90 * t85) * t85 - (t160 * t83 + t89 * t85) * t83;
t111 = (-t159 * t85 + t90 * t83) * t85 - (-t160 * t85 + t89 * t83) * t83;
t75 = rSges(2,1) * t85 - rSges(2,2) * t83;
t73 = rSges(2,1) * t83 + rSges(2,2) * t85;
t110 = -rSges(3,2) * t82 + t154;
t13 = -t104 * t83 + t40 * t85;
t14 = -t43 * t149 + t144;
t108 = t13 * t85 - t14 * t83;
t15 = -t40 * t83 + t113;
t35 = t45 * t146;
t16 = -t43 * t148 - t83 * t41 + t35;
t107 = t15 * t85 - t16 * t83;
t106 = -t23 * t83 + t24 * t85;
t76 = qJD(1) * t80;
t96 = t109 * qJD(2);
t31 = -qJD(1) * t122 + t76 + (-rSges(3,3) * qJD(1) - t96) * t85;
t32 = -t83 * t96 + (t110 * t85 - t153) * qJD(1);
t105 = -t31 * t85 - t32 * t83;
t102 = -t46 * t83 - t47 * t85;
t101 = -t84 * t68 - t82 * t70;
t92 = (t83 * t133 - t44 + t79) * t85 - (-t68 * t85 - t45) * t83;
t91 = (t140 * t82 - t141 * t84) * qJD(1);
t59 = t98 * qJD(2);
t60 = t99 * qJD(2);
t88 = -qJD(1) * t66 + t101 * qJD(2) + t59 * t82 - t60 * t84;
t87 = t158 * t84 + t92 * t82;
t86 = qJD(1) ^ 2;
t61 = t110 * qJD(2);
t54 = t109 * t83;
t22 = -t100 * t85 - t48;
t18 = t22 * qJD(1);
t17 = t102 * qJD(2);
t12 = -t61 * t127 - qJD(1) * t32 - qJDD(1) * t46 - t65 * t109 + (-qJDD(1) * t85 + t83 * t86) * pkin(1);
t11 = -t61 * t124 + qJD(1) * t31 + qJDD(1) * t47 + t64 * t109 + (-qJDD(1) * t83 - t85 * t86) * pkin(1);
t10 = -t161 * t85 + t88 * t83;
t9 = t161 * t83 + t88 * t85;
t8 = t103 * qJD(2) - t27 * t84 - t29 * t82;
t7 = t104 * qJD(2) - t28 * t84 - t30 * t82;
t6 = t107 * qJD(2) + t18;
t5 = t108 * qJD(2) + t126;
t1 = [(t18 + ((t113 - t14 + t144) * t85 + (-t13 - t35 + (-t104 + t41) * t83) * t83) * qJD(2)) * t120 - m(2) * (g(1) * t75 - g(3) * t73) + (t100 * qJD(2) + t59 * t84 + t60 * t82) * qJD(1) + (t20 + t22) * t157 + (t21 + t19) * t156 + (t8 + t9) * t118 + (t5 - t126 + ((t116 * t85 - t16 + t35) * t85 + (t116 * t83 + t144 - t15 - t38) * t83) * qJD(2)) * t117 + (t10 + t7 + t6) * t119 + (Icges(2,3) + t101 + m(2) * (t73 ^ 2 + t75 ^ 2)) * qJDD(1) + (-t23 * t76 + (t109 * t151 + t152) * qJD(2) + ((t154 * t23 + t155 * t24) * t83 + (-t110 * t24 + t155 * t23) * t85) * qJD(1) + (-g(3) + t12) * (-t155 * t85 + t115) + (-g(1) + t11) * (-t155 * t83 + t114)) * m(3); -t5 * t128 / 0.2e1 + t85 * (t10 * qJD(1) + t111 * qJD(2) + t21 * qJDD(1) + t13 * t65 + t14 * t64) / 0.2e1 + t108 * t156 + ((-t13 * t83 - t14 * t85) * qJD(1) + t111) * t119 + qJDD(1) * (t19 * t85 - t20 * t83) / 0.2e1 + qJD(1) * (t7 * t85 - t8 * t83 + (-t19 * t83 - t20 * t85) * qJD(1)) / 0.2e1 - t83 * (t9 * qJD(1) + t112 * qJD(2) + t22 * qJDD(1) + t15 * t65 + t16 * t64) / 0.2e1 + t107 * t157 + ((-t15 * t83 - t16 * t85) * qJD(1) + t112) * t118 + ((t127 * t48 - t125) * t85 + (t91 + (-t85 * t49 + t87) * qJD(2)) * t83) * t120 + ((t124 * t49 + t125) * t83 + (t91 + (-t83 * t48 + t87) * qJD(2)) * t85) * t117 - (t85 * t6 + (t140 * t84 + t141 * t82) * qJD(1) + (-t158 * t82 + t92 * t84) * qJD(2)) * qJD(1) / 0.2e1 + ((qJD(2) * t105 + t46 * t64 - t47 * t65) * t102 + t17 * ((-t46 * t85 + t47 * t83) * qJD(1) + t105) - t106 * t61 - (t11 * t83 + t12 * t85 + (-t23 * t85 - t151) * qJD(1)) * t109 - (t24 * t54 + t152) * qJD(1) - (t17 * (t54 * t83 + t55 * t85) - t106 * t110) * qJD(2) + g(1) * t54 + g(2) * t110 + g(3) * t55) * m(3);];
tau = t1;
