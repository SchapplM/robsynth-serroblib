% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:28
% EndTime: 2019-12-05 15:56:36
% DurationCPUTime: 2.70s
% Computational Cost: add. (2495->278), mult. (6715->409), div. (0->0), fcn. (5036->10), ass. (0->143)
t106 = cos(qJ(2));
t98 = sin(pkin(5));
t143 = t106 * t98;
t102 = sin(qJ(4));
t105 = cos(qJ(4));
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t83 = t102 * t97 - t105 * t99;
t111 = t83 * t143;
t151 = pkin(7) + qJ(3);
t88 = t151 * t97;
t89 = t151 * t99;
t55 = t102 * t89 + t105 * t88;
t180 = qJD(1) * t111 - qJD(3) * t83 - qJD(4) * t55;
t103 = sin(qJ(2));
t140 = qJD(1) * t98;
t135 = t103 * t140;
t78 = t83 * qJD(4);
t84 = t102 * t99 + t105 * t97;
t79 = t84 * qJD(4);
t186 = pkin(4) * t79 + pkin(8) * t78 - t135;
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t117 = Ifges(6,5) * t104 - Ifges(6,6) * t101;
t146 = Ifges(6,4) * t104;
t119 = -Ifges(6,2) * t101 + t146;
t147 = Ifges(6,4) * t101;
t121 = Ifges(6,1) * t104 - t147;
t122 = mrSges(6,1) * t101 + mrSges(6,2) * t104;
t145 = pkin(7) * qJD(2);
t87 = qJD(2) * qJ(3) + t135;
t100 = cos(pkin(5));
t138 = qJD(1) * t100;
t91 = t99 * t138;
t53 = t91 + (-t87 - t145) * t97;
t63 = t97 * t138 + t99 * t87;
t54 = t145 * t99 + t63;
t28 = t102 * t53 + t105 * t54;
t24 = qJD(4) * pkin(8) + t28;
t137 = qJD(1) * t106;
t133 = t98 * t137;
t124 = qJD(3) - t133;
t93 = -pkin(3) * t99 - pkin(2);
t71 = qJD(2) * t93 + t124;
t76 = t83 * qJD(2);
t77 = t84 * qJD(2);
t29 = pkin(4) * t76 - pkin(8) * t77 + t71;
t7 = -t101 * t24 + t104 * t29;
t8 = t101 * t29 + t104 * t24;
t127 = t8 * t101 + t7 * t104;
t163 = t104 / 0.2e1;
t164 = -t101 / 0.2e1;
t61 = qJD(4) * t101 + t104 * t77;
t167 = t61 / 0.2e1;
t161 = Ifges(6,4) * t61;
t179 = qJD(5) + t76;
t60 = qJD(4) * t104 - t101 * t77;
t17 = Ifges(6,2) * t60 + Ifges(6,6) * t179 + t161;
t59 = Ifges(6,4) * t60;
t18 = Ifges(6,1) * t61 + Ifges(6,5) * t179 + t59;
t27 = -t102 * t54 + t105 * t53;
t23 = -qJD(4) * pkin(4) - t27;
t185 = -t127 * mrSges(6,3) + t18 * t163 + t17 * t164 + t179 * t117 / 0.2e1 + t60 * t119 / 0.2e1 + t121 * t167 + t23 * t122;
t184 = t76 / 0.2e1;
t52 = pkin(4) * t83 - pkin(8) * t84 + t93;
t56 = -t102 * t88 + t105 * t89;
t22 = t101 * t52 + t104 * t56;
t183 = -qJD(5) * t22 - t180 * t101 + t186 * t104;
t21 = -t101 * t56 + t104 * t52;
t182 = qJD(5) * t21 + t186 * t101 + t180 * t104;
t150 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t60 + mrSges(6,2) * t61 + t77 * mrSges(5,3);
t112 = t84 * t143;
t181 = qJD(1) * t112 - qJD(3) * t84 - qJD(4) * t56;
t177 = t77 * Ifges(5,1) / 0.2e1;
t72 = Ifges(5,4) * t76;
t178 = t71 * mrSges(5,2) + Ifges(5,5) * qJD(4) - t72 / 0.2e1 + t177 + t185;
t176 = Ifges(5,2) * t184;
t148 = t97 ^ 2 + t99 ^ 2;
t175 = mrSges(4,3) * t148;
t82 = (qJD(3) + t133) * qJD(2);
t13 = qJD(4) * t27 - t82 * t83;
t136 = qJD(2) * t103;
t132 = t98 * t136;
t130 = qJD(1) * t132;
t69 = qJD(2) * t78;
t70 = qJD(2) * t79;
t35 = pkin(4) * t70 + pkin(8) * t69 + t130;
t1 = qJD(5) * t7 + t101 * t35 + t104 * t13;
t2 = -qJD(5) * t8 - t101 * t13 + t104 * t35;
t174 = t1 * t104 - t101 * t2;
t36 = qJD(5) * t60 - t104 * t69;
t37 = -qJD(5) * t61 + t101 * t69;
t173 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t36 + Ifges(6,6) * t37;
t172 = m(4) / 0.2e1;
t171 = t36 / 0.2e1;
t170 = t37 / 0.2e1;
t169 = -t60 / 0.2e1;
t168 = -t61 / 0.2e1;
t166 = t70 / 0.2e1;
t165 = -t179 / 0.2e1;
t162 = Ifges(5,4) * t77;
t14 = qJD(4) * t28 + t82 * t84;
t144 = t103 * t98;
t74 = t100 * t99 - t144 * t97;
t75 = t100 * t97 + t144 * t99;
t45 = t102 * t75 - t105 * t74;
t158 = t14 * t45;
t157 = t14 * t55;
t134 = t98 ^ 2 * t137;
t131 = t148 * t82;
t42 = t70 * mrSges(5,1) - t69 * mrSges(5,2);
t128 = -t1 * t101 - t2 * t104;
t126 = t7 * t101 - t8 * t104;
t125 = -(-t87 * t97 + t91) * t97 + t63 * t99;
t123 = mrSges(6,1) * t104 - mrSges(6,2) * t101;
t120 = Ifges(6,1) * t101 + t146;
t118 = Ifges(6,2) * t104 + t147;
t116 = Ifges(6,5) * t101 + Ifges(6,6) * t104;
t46 = t102 * t74 + t105 * t75;
t31 = -t101 * t46 - t104 * t143;
t115 = t101 * t143 - t104 * t46;
t113 = t125 * t106;
t110 = t7 * mrSges(6,1) + t71 * mrSges(5,1) + t179 * Ifges(6,3) + t61 * Ifges(6,5) + t60 * Ifges(6,6) - Ifges(5,6) * qJD(4) + t176 - t162 / 0.2e1 - t8 * mrSges(6,2);
t107 = qJD(2) ^ 2;
t86 = -qJD(2) * pkin(2) + t124;
t67 = Ifges(6,3) * t70;
t64 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t76;
t50 = pkin(4) * t77 + pkin(8) * t76;
t49 = mrSges(5,1) * t76 + mrSges(5,2) * t77;
t39 = mrSges(6,1) * t179 - mrSges(6,3) * t61;
t38 = -mrSges(6,2) * t179 + mrSges(6,3) * t60;
t26 = qJD(2) * t112 + qJD(4) * t46;
t25 = -qJD(2) * t111 - qJD(4) * t45;
t20 = -mrSges(6,2) * t70 + mrSges(6,3) * t37;
t19 = mrSges(6,1) * t70 - mrSges(6,3) * t36;
t15 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t12 = t101 * t50 + t104 * t27;
t11 = -t101 * t27 + t104 * t50;
t10 = t36 * Ifges(6,1) + t37 * Ifges(6,4) + t70 * Ifges(6,5);
t9 = t36 * Ifges(6,4) + t37 * Ifges(6,2) + t70 * Ifges(6,6);
t6 = qJD(5) * t115 - t101 * t25 + t104 * t132;
t5 = qJD(5) * t31 + t101 * t132 + t104 * t25;
t3 = [t45 * t15 + t31 * t19 - t115 * t20 + t25 * t64 + t5 * t38 + t6 * t39 + t150 * t26 + (-t45 * t69 - t46 * t70) * mrSges(5,3) + (-t106 * t42 + (t49 + qJD(2) * (-mrSges(4,1) * t99 + mrSges(4,2) * t97)) * t136) * t98 + m(5) * (t13 * t46 + t25 * t28 - t26 * t27 + t158) + m(6) * (-t1 * t115 + t2 * t31 + t23 * t26 + t5 * t8 + t6 * t7 + t158) + m(4) * (-t74 * t97 + t75 * t99) * t82 + 0.2e1 * (t98 * t113 * t172 + (m(5) * (t71 * t98 - t134) / 0.2e1 + (t86 * t98 - t134) * t172) * t103) * qJD(2) + (-mrSges(3,1) * t103 + (-mrSges(3,2) + t175) * t106) * t98 * t107; t55 * t15 + t21 * t19 + t22 * t20 + t93 * t42 + t180 * t64 + t183 * t39 + t182 * t38 - t49 * t135 + (t176 + t110) * t79 - (t177 + t178) * t78 + (t67 / 0.2e1 + mrSges(5,1) * t130 + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t70 + t173) * t83 + (t10 * t163 + t9 * t164 + t117 * t166 + t121 * t171 + t119 * t170 + t14 * t122 - Ifges(5,1) * t69 + mrSges(5,2) * t130 + t128 * mrSges(6,3) + (t23 * t123 + t116 * t165 + t118 * t169 + t120 * t168 - t104 * t17 / 0.2e1 + t18 * t164 + t126 * mrSges(6,3)) * qJD(5)) * t84 + (-t13 * t83 + t14 * t84 + t27 * t78 - t28 * t79 - t55 * t69 - t56 * t70) * mrSges(5,3) + (t124 * qJD(2) * t148 + t131) * mrSges(4,3) + (t83 * t69 - t70 * t84 + t78 * t184 - t77 * t79 / 0.2e1) * Ifges(5,4) - t150 * t181 + (t1 * t22 - t181 * t23 + t182 * t8 + t183 * t7 + t2 * t21 + t157) * m(6) + (t13 * t56 + t130 * t93 - t135 * t71 + t180 * t28 + t181 * t27 + t157) * m(5) + (-(t103 * t86 + t113) * t140 - pkin(2) * t130 + qJ(3) * t131 + t125 * qJD(3)) * m(4); t76 * t64 - t150 * t77 + (m(4) + m(5)) * t130 - t107 * t175 + (t179 * t38 + t19) * t104 + (-t179 * t39 + t20) * t101 - m(5) * (-t27 * t77 - t28 * t76) - m(4) * t125 * qJD(2) + t42 + (-t126 * t179 - t23 * t77 - t128) * m(6); t116 * t166 + t120 * t171 + t118 * t170 + t9 * t163 + t101 * t10 / 0.2e1 - Ifges(5,5) * t69 - Ifges(5,6) * t70 - t27 * t64 - t12 * t38 - t11 * t39 - t13 * mrSges(5,2) - pkin(4) * t15 - t150 * t28 + (-mrSges(5,1) - t123) * t14 + t174 * mrSges(6,3) + (t28 * mrSges(5,3) + t162 / 0.2e1 - t110) * t77 - (t72 / 0.2e1 + t27 * mrSges(5,3) + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t77 - t178) * t76 + t185 * qJD(5) + (-pkin(4) * t14 - t11 * t7 - t12 * t8 - t23 * t28) * m(6) + (-t101 * t19 + t104 * t20 + m(6) * t174 + (-m(6) * t127 - t101 * t38 - t104 * t39) * qJD(5)) * pkin(8); t67 - t23 * (mrSges(6,1) * t61 + mrSges(6,2) * t60) + (Ifges(6,1) * t60 - t161) * t168 + t17 * t167 + (Ifges(6,5) * t60 - Ifges(6,6) * t61) * t165 - t7 * t38 + t8 * t39 + (t60 * t7 + t61 * t8) * mrSges(6,3) + (-Ifges(6,2) * t61 + t18 + t59) * t169 + t173;];
tauc = t3(:);
