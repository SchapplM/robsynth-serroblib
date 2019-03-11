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
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
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
% StartTime: 2019-03-08 18:00:00
% EndTime: 2019-03-08 18:00:02
% DurationCPUTime: 1.68s
% Computational Cost: add. (1143->211), mult. (3010->317), div. (0->0), fcn. (2394->4), ass. (0->131)
t87 = sin(qJ(2));
t89 = cos(qJ(2));
t101 = Icges(3,5) * t87 + Icges(3,6) * t89;
t88 = sin(qJ(1));
t49 = t101 * t88;
t90 = cos(qJ(1));
t50 = t101 * t90;
t117 = rSges(3,1) * t87 + rSges(3,2) * t89;
t132 = t90 * qJD(2);
t176 = t117 * t132;
t135 = qJD(2) * t88;
t175 = t117 * t135;
t146 = Icges(3,4) * t89;
t104 = -Icges(3,2) * t87 + t146;
t42 = Icges(3,6) * t90 + t104 * t88;
t91 = qJD(1) ^ 2;
t143 = Icges(3,5) * t90;
t148 = Icges(3,1) * t89;
t147 = Icges(3,4) * t87;
t80 = t88 * t147;
t44 = t148 * t88 + t143 - t80;
t162 = t44 * t89;
t112 = t42 * t87 - t162;
t138 = Icges(3,3) * t90;
t140 = Icges(3,6) * t88;
t144 = Icges(3,5) * t89;
t40 = -t140 * t87 + t144 * t88 + t138;
t13 = -t112 * t88 + t90 * t40;
t102 = -Icges(3,6) * t87 + t144;
t103 = Icges(3,2) * t89 + t147;
t105 = Icges(3,1) * t87 + t146;
t107 = -t103 * t87 + t105 * t89;
t173 = t107 * qJD(1) + t102 * qJD(2);
t145 = Icges(3,5) * t88;
t81 = t90 * t147;
t45 = -t148 * t90 + t145 + t81;
t161 = t45 * t89;
t43 = -t104 * t90 + t140;
t110 = t43 * t87 - t161;
t98 = qJD(2) * t101;
t172 = -t90 * t98 + (-t102 * t88 + t110 - t138) * qJD(1);
t41 = Icges(3,3) * t88 - t102 * t90;
t137 = qJD(1) * t41;
t171 = qJD(1) * t112 + t88 * t98 + t137;
t170 = (t45 + t81) * t88 + (t44 - t80) * t90;
t131 = qJD(1) * qJD(2);
t63 = qJDD(2) * t88 + t131 * t90;
t169 = t63 / 0.2e1;
t64 = -qJDD(2) * t90 + t131 * t88;
t168 = t64 / 0.2e1;
t167 = rSges(3,3) + pkin(1);
t165 = rSges(3,2) * t87;
t158 = t89 * t88;
t130 = rSges(3,1) * t158;
t160 = t87 * t88;
t82 = rSges(3,2) * t160;
t122 = -t82 + t130;
t46 = rSges(3,3) * t90 + t122;
t136 = qJD(1) * t90;
t85 = pkin(1) * t136;
t24 = qJD(1) * t46 + t176 + t85;
t163 = t24 * t90;
t159 = t87 * t90;
t86 = t88 * rSges(3,3);
t157 = t89 * t90;
t154 = t42 * t159 + t88 * t40;
t153 = t43 * t159 + t88 * t41;
t106 = -t147 + t148;
t150 = -t103 + t106;
t149 = -t104 - t105;
t21 = t103 * t160 - t105 * t158 - t50;
t134 = t21 * qJD(1);
t133 = t102 * qJD(1);
t47 = -rSges(3,1) * t157 + rSges(3,2) * t159 + t86;
t36 = pkin(1) * t88 + t47;
t129 = -t135 / 0.2e1;
t128 = t135 / 0.2e1;
t127 = -t132 / 0.2e1;
t126 = t132 / 0.2e1;
t33 = t43 * t160;
t125 = t90 * t41 - t33;
t124 = t40 + t161;
t123 = rSges(3,3) * t136 + qJD(1) * t130 + t176;
t111 = t43 * t89 + t45 * t87;
t99 = qJD(2) * t103;
t27 = qJD(1) * t42 + t90 * t99;
t100 = qJD(2) * t105;
t29 = t90 * t100 + (t106 * t88 + t143) * qJD(1);
t94 = qJD(2) * t111 + t27 * t87 - t29 * t89 + t137;
t19 = t42 * t89 + t44 * t87;
t28 = qJD(1) * t43 + t88 * t99;
t30 = t88 * t100 + (-t106 * t90 + t145) * qJD(1);
t95 = -qJD(1) * t40 - qJD(2) * t19 + t28 * t87 - t30 * t89;
t121 = -(t171 * t88 + t95 * t90) * t90 + (-t172 * t88 + t94 * t90) * t88;
t120 = -(-t171 * t90 + t95 * t88) * t90 + (t172 * t90 + t94 * t88) * t88;
t119 = rSges(2,1) * t90 - rSges(2,2) * t88;
t72 = rSges(2,1) * t88 + rSges(2,2) * t90;
t118 = rSges(3,1) * t89 - t165;
t14 = -t158 * t45 - t125;
t116 = -t13 * t90 + t14 * t88;
t15 = t157 * t44 - t154;
t16 = -t157 * t45 + t153;
t115 = -t15 * t90 + t16 * t88;
t23 = qJD(1) * t36 + t175;
t114 = -t23 * t88 - t163;
t31 = -qJD(1) * t82 + t123;
t55 = t117 * t88;
t32 = qJD(2) * t55 + (-t118 * t90 + t86) * qJD(1);
t113 = t31 * t90 + t32 * t88;
t109 = -t46 * t88 + t47 * t90;
t108 = -t103 * t89 - t105 * t87;
t97 = t42 * t90 + t43 * t88;
t96 = (t149 * t87 + t150 * t89) * qJD(1);
t58 = t104 * qJD(2);
t59 = t106 * qJD(2);
t93 = -qJD(1) * t101 + qJD(2) * t108 - t58 * t87 + t59 * t89;
t92 = t170 * t87 + t97 * t89;
t60 = t118 * qJD(2);
t56 = t117 * t90;
t22 = t107 * t90 - t49;
t18 = t22 * qJD(1);
t17 = t109 * qJD(2);
t12 = t60 * t132 - qJD(1) * t32 + qJDD(1) * t46 - t64 * t117 + (qJDD(1) * t90 - t88 * t91) * pkin(1);
t11 = t60 * t135 + qJD(1) * t31 + qJDD(1) * t47 + t63 * t117 + (qJDD(1) * t88 + t90 * t91) * pkin(1);
t10 = t173 * t90 + t93 * t88;
t9 = -t173 * t88 + t93 * t90;
t8 = qJD(2) * t110 - t27 * t89 - t29 * t87;
t7 = -qJD(2) * t112 - t28 * t89 - t30 * t87;
t6 = qJD(2) * t115 + t18;
t5 = qJD(2) * t116 - t134;
t1 = [(t18 + ((-t33 + t14 + (t41 - t162) * t90 + t154) * t90 + t153 * t88) * qJD(2)) * t126 - t64 * t21 / 0.2e1 + t19 * t168 + (t107 * qJD(2) + t58 * t89 + t59 * t87) * qJD(1) + (-t111 + t22) * t169 + (t5 + t134 + ((t124 * t90 - t153 + t16) * t90 + (t124 * t88 + t125 + t15) * t88) * qJD(2)) * t129 + (t8 + t9) * t128 + (Icges(2,3) - t108) * qJDD(1) + (t23 * (t85 + t123) - t24 * t175 + (t118 * t163 + (-t165 * t23 - t167 * t24) * t88) * qJD(1) + (-g(3) + t12) * (t167 * t90 + t122) + (-g(1) + t11) * t36) * m(3) + ((qJDD(1) * t119 + g(1)) * t119 + (qJDD(1) * t72 - g(3)) * t72) * m(2) + (t10 + t7 + t6) * t127; -t90 * (t10 * qJD(1) + t120 * qJD(2) - t21 * qJDD(1) + t13 * t64 + t14 * t63) / 0.2e1 + t116 * t168 + ((t13 * t88 + t14 * t90) * qJD(1) + t120) * t127 + qJDD(1) * (-t111 * t88 - t19 * t90) / 0.2e1 + t6 * t136 / 0.2e1 + t88 * (t9 * qJD(1) + t121 * qJD(2) + t22 * qJDD(1) + t15 * t64 + t16 * t63) / 0.2e1 + t115 * t169 + ((t15 * t88 + t16 * t90) * qJD(1) + t121) * t128 + ((t132 * t49 + t133) * t90 + (t96 + (-t90 * t50 + t92) * qJD(2)) * t88) * t126 - qJD(1) * ((-t149 * t89 + t150 * t87) * qJD(1) + (-t170 * t89 + t97 * t87) * qJD(2)) / 0.2e1 + ((t135 * t50 - t133) * t88 + (t96 + (-t88 * t49 + t92) * qJD(2)) * t90) * t129 + ((-qJD(1) * t111 - t7) * t90 + (qJD(1) * t19 + t5 + t8) * t88) * qJD(1) / 0.2e1 + ((qJD(2) * t113 - t46 * t63 - t47 * t64) * t109 + t17 * ((-t46 * t90 - t47 * t88) * qJD(1) + t113) - t114 * t60 - (-t11 * t88 - t12 * t90 + (-t23 * t90 + t24 * t88) * qJD(1)) * t117 - (t23 * t56 - t24 * t55) * qJD(1) - (t17 * (t55 * t88 + t56 * t90) - t114 * t118) * qJD(2) - g(1) * t55 + g(2) * t118 - g(3) * t56) * m(3);];
tau  = t1;
