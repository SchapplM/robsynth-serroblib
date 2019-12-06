% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:47
% EndTime: 2019-12-05 16:50:56
% DurationCPUTime: 2.78s
% Computational Cost: add. (1841->271), mult. (4707->371), div. (0->0), fcn. (2955->6), ass. (0->138)
t180 = Ifges(5,1) + Ifges(6,1);
t185 = Ifges(6,4) + Ifges(5,5);
t186 = mrSges(5,1) + mrSges(6,1);
t184 = -Ifges(5,6) + Ifges(6,6);
t106 = qJD(3) + qJD(4);
t109 = sin(qJ(4));
t112 = cos(qJ(4));
t113 = cos(qJ(3));
t141 = qJD(2) * t113;
t110 = sin(qJ(3));
t143 = qJD(2) * t110;
t77 = t109 * t143 - t112 * t141;
t160 = t77 * Ifges(6,5);
t74 = Ifges(5,4) * t77;
t146 = t109 * t113;
t84 = t110 * t112 + t146;
t78 = t84 * qJD(2);
t183 = t185 * t106 + t180 * t78 + t160 - t74;
t111 = sin(qJ(2));
t175 = t111 * t106;
t182 = -Ifges(4,1) / 0.2e1;
t181 = -Ifges(4,4) * t141 / 0.2e1;
t157 = mrSges(6,2) + mrSges(5,3);
t156 = -Ifges(5,4) + Ifges(6,5);
t114 = cos(qJ(2));
t147 = t109 * t110;
t83 = -t112 * t113 + t147;
t118 = t83 * t114;
t169 = -pkin(7) - pkin(6);
t131 = qJD(3) * t169;
t123 = t113 * t131;
t92 = t169 * t110;
t93 = t169 * t113;
t57 = -t109 * t93 - t112 * t92;
t87 = t110 * t131;
t179 = qJD(1) * t118 - qJD(4) * t57 + t109 * t123 + t112 * t87;
t144 = qJD(1) * t114;
t58 = t109 * t92 - t112 * t93;
t178 = qJD(4) * t58 + t109 * t87 - t112 * t123 - t84 * t144;
t47 = mrSges(6,1) * t77 - mrSges(6,3) * t78;
t48 = mrSges(5,1) * t77 + mrSges(5,2) * t78;
t177 = -t47 - t48;
t145 = qJD(1) * t111;
t95 = qJD(2) * pkin(6) + t145;
t127 = pkin(7) * qJD(2) + t95;
t70 = t127 * t113;
t152 = t109 * t70;
t69 = t127 * t110;
t66 = qJD(3) * pkin(3) - t69;
t27 = t112 * t66 - t152;
t176 = qJD(5) - t27;
t107 = t110 ^ 2;
t108 = t113 ^ 2;
t173 = -t77 / 0.2e1;
t172 = t77 / 0.2e1;
t170 = t78 / 0.2e1;
t138 = qJD(3) * t113;
t133 = t95 * t138;
t117 = -t133 + (-pkin(7) * t138 - t110 * t144) * qJD(2);
t150 = t112 * t70;
t28 = t109 * t66 + t150;
t136 = qJD(1) * qJD(2);
t128 = t114 * t136;
t94 = t113 * t128;
t55 = -qJD(3) * t69 + t94;
t5 = qJD(4) * t28 + t109 * t55 - t112 * t117;
t168 = t5 * t57;
t71 = t84 * t111;
t167 = t5 * t71;
t166 = t5 * t84;
t165 = -t106 / 0.2e1;
t163 = mrSges(6,2) * t78;
t162 = mrSges(5,3) * t77;
t20 = -pkin(4) * t106 + t176;
t161 = t20 * t77;
t159 = t78 * mrSges(5,3);
t158 = t78 * Ifges(5,4);
t59 = -mrSges(5,2) * t106 - t162;
t62 = -t77 * mrSges(6,2) + mrSges(6,3) * t106;
t155 = t59 + t62;
t154 = t186 * t106 - t159 - t163;
t153 = Ifges(4,4) * t110;
t96 = -qJD(2) * pkin(2) - t144;
t151 = t111 * t96;
t101 = t111 * t136;
t139 = qJD(3) * t110;
t134 = pkin(3) * t139;
t82 = qJD(2) * t134 + t101;
t149 = Ifges(4,5) * qJD(3);
t148 = Ifges(4,6) * qJD(3);
t142 = qJD(2) * t111;
t140 = qJD(2) * t114;
t137 = qJD(4) * t112;
t135 = qJD(2) * qJD(3);
t104 = -pkin(3) * t113 - pkin(2);
t130 = t149 / 0.2e1;
t129 = -t148 / 0.2e1;
t126 = (t107 + t108) * t95;
t125 = t148 / 0.2e1 + (t113 * Ifges(4,2) + t153) * qJD(2) / 0.2e1 - t96 * mrSges(4,1);
t124 = t143 * t182 + t181 - t149 / 0.2e1 - t96 * mrSges(4,2);
t44 = pkin(4) * t78 + qJ(5) * t77;
t63 = -t139 * t95 + t94;
t64 = -t110 * t128 - t133;
t121 = -t64 * t110 + t63 * t113;
t90 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t143;
t91 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t141;
t120 = t110 * t91 + t113 * t90;
t119 = t110 * t90 - t113 * t91;
t4 = -qJD(4) * t152 + t109 * t117 + t112 * t55 + t66 * t137;
t81 = qJD(2) * t104 - t144;
t54 = t106 * t84;
t53 = t106 * t83;
t2 = qJD(5) * t106 + t4;
t25 = qJ(5) * t106 + t28;
t26 = pkin(4) * t77 - qJ(5) * t78 + t81;
t73 = Ifges(6,5) * t78;
t34 = t106 * Ifges(6,6) + t77 * Ifges(6,3) + t73;
t35 = -t77 * Ifges(5,2) + t106 * Ifges(5,6) + t158;
t42 = t53 * qJD(2);
t43 = t54 * qJD(2);
t116 = -t4 * mrSges(5,2) + t2 * mrSges(6,3) - t27 * t162 + t25 * t163 - t26 * (mrSges(6,1) * t78 + mrSges(6,3) * t77) - t81 * (mrSges(5,1) * t78 - mrSges(5,2) * t77) + t35 * t170 + (Ifges(6,3) * t78 - t160) * t173 - t186 * t5 + t184 * t43 - t185 * t42 + (t184 * t78 - t185 * t77) * t165 + (-Ifges(5,2) * t78 + t183 - t74) * t172 - (-t180 * t77 - t158 + t34 + t73) * t78 / 0.2e1;
t115 = qJD(2) ^ 2;
t103 = -pkin(3) * t112 - pkin(4);
t100 = pkin(3) * t109 + qJ(5);
t97 = pkin(3) * t137 + qJD(5);
t80 = (mrSges(4,1) * t110 + mrSges(4,2) * t113) * t135;
t72 = t83 * t111;
t49 = pkin(4) * t83 - qJ(5) * t84 + t104;
t31 = pkin(3) * t143 + t44;
t30 = -t112 * t69 - t152;
t29 = -t109 * t69 + t150;
t13 = t140 * t146 + (t110 * t140 + t113 * t175) * t112 - t147 * t175;
t12 = -qJD(2) * t118 - t84 * t175;
t9 = mrSges(5,1) * t43 - mrSges(5,2) * t42;
t8 = mrSges(6,1) * t43 + mrSges(6,3) * t42;
t7 = pkin(4) * t54 + qJ(5) * t53 - qJD(5) * t84 + t134;
t6 = pkin(4) * t43 + qJ(5) * t42 - qJD(5) * t78 + t82;
t1 = [-t154 * t13 + t155 * t12 + (-t115 * mrSges(3,2) - t119 * qJD(2) - t8 - t80 - t9) * t114 + (-t115 * mrSges(3,1) - t120 * qJD(3) + (qJD(2) * (-mrSges(4,1) * t113 + mrSges(4,2) * t110) - t177) * qJD(2)) * t111 + m(4) * (t121 * t111 + (t151 + (t126 - t145) * t114) * qJD(2)) + m(5) * (-t114 * t82 + t12 * t28 - t13 * t27 + t142 * t81 - t4 * t72 + t167) + m(6) * (-t114 * t6 + t12 * t25 + t13 * t20 + t142 * t26 - t2 * t72 + t167) + t157 * (-t42 * t71 + t43 * t72); (0.3e1 / 0.2e1 * t108 - 0.3e1 / 0.2e1 * t107) * Ifges(4,4) * t135 + (-t2 * t83 + t166) * mrSges(6,2) + (-t4 * t83 + t166) * mrSges(5,3) - (t156 * t83 + t157 * t57 + t180 * t84) * t42 + (t156 * t84 + (Ifges(5,2) + Ifges(6,3)) * t83 - t157 * t58) * t43 + ((-pkin(6) * t90 - t124 + t130) * t113 + (-pkin(6) * t91 + pkin(3) * t48 + t129 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t141 - t125) * t110) * qJD(3) + (t177 * t111 + t119 * t114) * qJD(1) + t121 * mrSges(4,3) + t104 * t9 - pkin(2) * t80 + t6 * (mrSges(6,1) * t83 - mrSges(6,3) * t84) + t82 * (mrSges(5,1) * t83 + mrSges(5,2) * t84) + t7 * t47 + t49 * t8 - (t81 * mrSges(5,2) + t20 * mrSges(6,2) - t27 * mrSges(5,3) - t26 * mrSges(6,3) + Ifges(5,4) * t173 + Ifges(6,5) * t172) * t53 + (-t25 * mrSges(6,2) - t28 * mrSges(5,3) + Ifges(6,3) * t172 - Ifges(5,2) * t173 + t81 * mrSges(5,1) + t26 * mrSges(6,1) + t34 / 0.2e1 - t35 / 0.2e1) * t54 - t183 * t53 / 0.2e1 + (t156 * t54 - t180 * t53) * t170 + (t184 * t54 - t185 * t53) * t106 / 0.2e1 + t179 * t155 - t178 * t154 + (t2 * t58 + t49 * t6 + t168 + (-t145 + t7) * t26 + t179 * t25 + t178 * t20) * m(6) + (t104 * t82 + t4 * t58 + t168 + (t134 - t145) * t81 + t179 * t28 - t178 * t27) * m(5) + (-(t114 * t126 + t151) * qJD(1) - pkin(2) * t101 + t121 * pkin(6)) * m(4); ((t130 + t181 + t124) * t113 + (t129 + (t153 / 0.2e1 + (t182 + Ifges(4,2) / 0.2e1) * t113) * qJD(2) + (-m(5) * t81 - t48) * pkin(3) + t125) * t110) * qJD(2) + (m(5) * (t109 * t4 - t112 * t5) + (-t109 * t43 + t112 * t42) * mrSges(5,3) + ((m(5) * t28 + t59) * t112 + (-m(5) * t27 + m(6) * t20 - t154) * t109) * qJD(4)) * pkin(3) - m(5) * (-t27 * t29 + t28 * t30) + (-t100 * t43 - t103 * t42 + t161) * mrSges(6,2) + t97 * t62 - t63 * mrSges(4,2) + t64 * mrSges(4,1) - t31 * t47 + t28 * t159 - m(6) * (t20 * t29 + t25 * t30 + t26 * t31) + t116 + t120 * t95 - t155 * t30 + t154 * t29 + m(6) * (t100 * t2 + t103 * t5 + t25 * t97); qJD(5) * t62 - t44 * t47 + t116 + (t154 + t159) * t28 - t155 * t27 + (pkin(4) * t42 - qJ(5) * t43 + t161) * mrSges(6,2) + (-pkin(4) * t5 + qJ(5) * t2 + t176 * t25 - t20 * t28 - t26 * t44) * m(6); -t42 * mrSges(6,2) - t106 * t62 + t78 * t47 + 0.2e1 * (t5 / 0.2e1 + t25 * t165 + t26 * t170) * m(6);];
tauc = t1(:);
