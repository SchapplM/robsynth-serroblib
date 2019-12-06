% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:38
% EndTime: 2019-12-05 18:08:49
% DurationCPUTime: 3.22s
% Computational Cost: add. (1702->354), mult. (4773->498), div. (0->0), fcn. (3118->6), ass. (0->166)
t204 = Ifges(5,2) / 0.2e1;
t126 = qJD(3) * qJD(4);
t77 = cos(qJ(4));
t130 = qJD(4) * t77;
t78 = cos(qJ(3));
t133 = qJD(3) * t78;
t74 = sin(qJ(4));
t75 = sin(qJ(3));
t42 = t74 * t126 + (t75 * t130 + t74 * t133) * qJD(1);
t181 = t42 / 0.2e1;
t73 = sin(qJ(5));
t176 = -t73 / 0.2e1;
t137 = qJD(1) * t78;
t67 = qJD(4) - t137;
t203 = Ifges(5,6) * t67 / 0.2e1;
t202 = qJD(3) / 0.2e1;
t138 = qJD(1) * t75;
t62 = qJD(3) * t74 + t77 * t138;
t170 = Ifges(5,4) * t62;
t134 = qJD(3) * t77;
t61 = -t74 * t138 + t134;
t201 = t203 + t61 * t204 + t170 / 0.2e1;
t120 = qJD(3) * t138;
t200 = t120 / 0.2e1;
t199 = t138 / 0.2e1;
t118 = -Ifges(4,6) * qJD(3) / 0.2e1;
t119 = Ifges(4,5) * t202;
t166 = Ifges(5,5) * t67;
t129 = qJ(2) * qJD(1);
t121 = t78 * t129;
t59 = -t77 * qJD(2) + t74 * t121;
t76 = cos(qJ(5));
t108 = mrSges(6,1) * t73 + mrSges(6,2) * t76;
t94 = mrSges(5,3) + t108;
t198 = t94 * t59;
t122 = t75 * t129;
t197 = mrSges(5,1) * t122;
t196 = mrSges(5,2) * t122;
t132 = qJD(4) * t59;
t128 = qJ(2) * qJD(3);
t57 = (qJD(2) * t78 - t75 * t128) * qJD(1);
t28 = t57 * t77 - t132;
t131 = qJD(4) * t74;
t125 = t75 * t131;
t41 = t77 * t126 + (t77 * t133 - t125) * qJD(1);
t195 = t28 * mrSges(5,2) - Ifges(5,5) * t41 + Ifges(5,6) * t42;
t154 = t59 * t74;
t135 = qJD(3) * mrSges(4,2);
t64 = mrSges(4,3) * t137 - t135;
t60 = qJD(2) * t74 + t77 * t121;
t99 = t60 * t77 + t154;
t89 = -m(5) * t99 - m(6) * t154 - t64;
t35 = t62 * t76 + t67 * t73;
t17 = -t35 * qJD(5) + t76 * t120 - t41 * t73;
t14 = Ifges(6,6) * t17;
t34 = -t62 * t73 + t67 * t76;
t16 = t34 * qJD(5) + t73 * t120 + t41 * t76;
t15 = Ifges(6,5) * t16;
t1 = Ifges(6,3) * t42 + t14 + t15;
t36 = t76 * t122 - t60 * t73;
t58 = (qJD(2) * t75 + t78 * t128) * qJD(1);
t7 = t36 * qJD(5) + t28 * t76 + t58 * t73;
t37 = t73 * t122 + t60 * t76;
t8 = -t37 * qJD(5) - t28 * t73 + t58 * t76;
t111 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t194 = (Ifges(6,3) / 0.2e1 + t204) * t42 - t28 * mrSges(5,3) + t1 / 0.2e1 + Ifges(5,2) * t181 - Ifges(5,6) * t200 - t41 * Ifges(5,4) + t15 / 0.2e1 + t14 / 0.2e1 + t58 * mrSges(5,1) + t111;
t193 = -t36 * mrSges(6,1) + t37 * mrSges(6,2) + t201;
t114 = -qJD(5) + t134;
t192 = t114 * t75 + t78 * t131;
t136 = qJD(3) * mrSges(4,1);
t142 = mrSges(5,1) * t61 - mrSges(5,2) * t62 - mrSges(4,3) * t138 + t136;
t69 = Ifges(4,4) * t137;
t191 = -t142 * qJ(2) + qJD(2) * mrSges(4,2) + Ifges(4,1) * t199 + t119 + t69 / 0.2e1;
t169 = Ifges(6,4) * t35;
t56 = qJD(5) - t61;
t10 = Ifges(6,2) * t34 + Ifges(6,6) * t56 + t169;
t101 = t36 * t76 + t37 * t73;
t190 = -t101 * mrSges(6,3) + t10 * t176;
t3 = Ifges(6,1) * t16 + Ifges(6,4) * t17 + Ifges(6,5) * t42;
t189 = t3 / 0.2e1;
t188 = -t10 / 0.2e1;
t187 = t16 / 0.2e1;
t186 = t17 / 0.2e1;
t185 = -t34 / 0.2e1;
t184 = t34 / 0.2e1;
t183 = -t35 / 0.2e1;
t182 = t35 / 0.2e1;
t180 = -t56 / 0.2e1;
t179 = t56 / 0.2e1;
t175 = t76 / 0.2e1;
t174 = m(5) * qJ(2) ^ 2;
t173 = -mrSges(5,1) * t120 - mrSges(6,1) * t17 + mrSges(6,2) * t16 + mrSges(5,3) * t41;
t172 = mrSges(5,3) * t62;
t171 = Ifges(5,1) * t62;
t55 = Ifges(5,4) * t61;
t168 = Ifges(6,4) * t73;
t167 = Ifges(6,4) * t76;
t165 = Ifges(6,5) * t35;
t162 = Ifges(6,6) * t34;
t161 = Ifges(6,3) * t56;
t149 = t75 * t76;
t148 = t75 * t77;
t33 = Ifges(6,4) * t34;
t11 = Ifges(6,1) * t35 + Ifges(6,5) * t56 + t33;
t147 = t76 * t11;
t146 = t76 * t78;
t43 = -mrSges(5,2) * t67 + t61 * mrSges(5,3);
t145 = t77 * t43;
t144 = t77 * t78;
t143 = -mrSges(5,1) * t67 - mrSges(6,1) * t34 + mrSges(6,2) * t35 + t172;
t79 = qJD(1) ^ 2;
t141 = qJ(2) * t79;
t127 = qJD(1) * qJD(2);
t117 = m(3) * qJ(2) + mrSges(3,3);
t116 = qJ(2) * t127;
t115 = -qJD(5) * t77 + qJD(3);
t110 = t7 * t76 - t8 * t73;
t109 = mrSges(6,1) * t76 - mrSges(6,2) * t73;
t107 = Ifges(6,1) * t76 - t168;
t106 = Ifges(6,1) * t73 + t167;
t105 = -Ifges(6,2) * t73 + t167;
t104 = Ifges(6,2) * t76 + t168;
t103 = Ifges(6,5) * t76 - Ifges(6,6) * t73;
t102 = Ifges(6,5) * t73 + Ifges(6,6) * t76;
t100 = t36 * t73 - t37 * t76;
t98 = t115 * t78;
t97 = t76 * t144 + t73 * t75;
t54 = t76 * t148 - t73 * t78;
t96 = -t73 * t144 + t149;
t95 = t73 * t148 + t146;
t23 = -mrSges(6,2) * t56 + mrSges(6,3) * t34;
t24 = mrSges(6,1) * t56 - mrSges(6,3) * t35;
t93 = t23 * t76 - t24 * t73 + t43;
t29 = t60 * qJD(4) + t57 * t74;
t92 = t59 * t130 + t29 * t74;
t91 = t58 * mrSges(5,2) + t41 * Ifges(5,1) - t42 * Ifges(5,4) + Ifges(5,5) * t200;
t27 = t166 + t55 + t171;
t90 = t27 / 0.2e1 + t171 / 0.2e1 + t55 / 0.2e1 + t166 / 0.2e1;
t88 = t59 * mrSges(5,3) + t90;
t87 = -t88 - t196;
t86 = qJD(2) * mrSges(4,1) + t67 * Ifges(5,3) + t62 * Ifges(5,5) + t61 * Ifges(5,6) + t118 - (Ifges(4,4) * t75 + t78 * Ifges(4,2)) * qJD(1) / 0.2e1 - t59 * mrSges(5,1) - t60 * mrSges(5,2);
t9 = t161 + t162 + t165;
t84 = t161 / 0.2e1 + t162 / 0.2e1 + t165 / 0.2e1 + t9 / 0.2e1 - t60 * mrSges(5,3) - t193 - t201;
t83 = t84 * t74;
t82 = -t84 - t197;
t81 = t105 * t184 + t107 * t182 + t103 * t179 + t147 / 0.2e1 + t190;
t72 = t78 ^ 2;
t71 = t75 ^ 2;
t68 = t71 * t141;
t66 = Ifges(5,3) * t120;
t65 = t71 * t116;
t50 = t97 * qJ(2);
t49 = t96 * qJ(2);
t48 = t97 * qJD(1);
t47 = t96 * qJD(1);
t46 = t54 * t129;
t45 = t95 * t129;
t31 = -mrSges(5,2) * t120 - mrSges(5,3) * t42;
t22 = t114 * t146 + (t115 * t73 - t76 * t131) * t75;
t21 = t115 * t149 + (-t114 * t78 + t125) * t73;
t13 = t96 * qJD(2) + (t192 * t73 + t76 * t98) * qJ(2);
t12 = t97 * qJD(2) + (-t192 * t76 + t73 * t98) * qJ(2);
t6 = -mrSges(6,2) * t42 + mrSges(6,3) * t17;
t5 = mrSges(6,1) * t42 - mrSges(6,3) * t16;
t2 = Ifges(6,4) * t16 + Ifges(6,2) * t17 + Ifges(6,6) * t42;
t4 = [(Ifges(6,1) * t22 + Ifges(6,4) * t21) * t182 + (Ifges(6,1) * t54 - Ifges(6,4) * t95) * t187 + t59 * (-mrSges(6,1) * t21 + mrSges(6,2) * t22) - t95 * t2 / 0.2e1 + t29 * (mrSges(6,1) * t95 + mrSges(6,2) * t54) + t54 * t189 + t49 * t5 + t50 * t6 + t21 * t10 / 0.2e1 + t22 * t11 / 0.2e1 + t12 * t23 + t13 * t24 + (Ifges(6,5) * t22 + Ifges(6,6) * t21) * t179 + (Ifges(6,5) * t54 - Ifges(6,6) * t95) * t181 + (Ifges(6,4) * t22 + Ifges(6,2) * t21) * t184 + (Ifges(6,4) * t54 - Ifges(6,2) * t95) * t186 + 0.2e1 * t117 * t127 + m(6) * (t12 * t37 + t13 * t36 + t49 * t8 + t50 * t7) + m(5) * t65 + m(4) * (t72 * t116 + t65) + (t21 * t37 - t22 * t36 - t54 * t8 - t7 * t95) * mrSges(6,3) + (t57 * mrSges(4,3) - t66 / 0.2e1 + t29 * mrSges(5,1) + (t143 * t74 + t145 - t89) * qJD(2) + (t77 * t31 + t173 * t74 + (t143 * t77 - t74 * t43) * qJD(4) + m(6) * t92 + m(4) * t57 + m(5) * (-t60 * t131 + t28 * t77 + t92)) * qJ(2) + (0.3e1 / 0.2e1 * t69 + t119 + t88 * t77 + t83 + t191) * qJD(3) + t195) * t78 + (-t142 * qJD(2) + qJ(2) * (mrSges(5,1) * t42 + mrSges(5,2) * t41) + (t29 * mrSges(5,3) + t91) * t77 + t194 * t74 + (t87 * t74 - t82 * t77) * qJD(4) + (t118 + ((Ifges(5,5) * t77 / 0.2e1 - Ifges(5,6) * t74 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t75 + (t174 - Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t78) * qJD(1) + ((mrSges(5,2) * t137 - t43) * t77 + (mrSges(5,1) * t137 - t143) * t74 + t89) * qJ(2) + t86) * qJD(3) + (mrSges(4,3) + (m(5) + m(4)) * qJ(2)) * t58) * t75; -t48 * t23 - t47 * t24 - t117 * t79 + (t93 * qJD(4) - t173) * t77 + ((t135 - t64 - t145) * t78 + (t136 + t142) * t75) * qJD(1) + (-t73 * t5 + t76 * t6 + t31 + (-t23 * t73 - t24 * t76) * qJD(5) + t67 * t143) * t74 - m(4) * (t72 * t141 + t68) + (-t137 * t154 - t36 * t47 - t37 * t48 + (-t100 * qJD(4) - t29) * t77 + (-t101 * qJD(5) + t110 + t132) * t74) * m(6) + (t28 * t74 - t29 * t77 + t67 * t99 - t68) * m(5); (Ifges(6,5) * t48 + Ifges(6,6) * t47) * t180 + (Ifges(6,4) * t48 + Ifges(6,2) * t47) * t185 + (Ifges(6,1) * t48 + Ifges(6,4) * t47) * t183 - t58 * mrSges(4,1) - t59 * (-mrSges(6,1) * t47 + mrSges(6,2) * t48) - t57 * mrSges(4,2) - t45 * t24 + t46 * t23 + t47 * t188 - t48 * t11 / 0.2e1 + (t36 * t48 - t37 * t47) * mrSges(6,3) - m(6) * (t36 * t45 - t37 * t46) - t75 * t78 * t79 * t174 - t194 * t77 + (t107 * t187 + t105 * t186 + t103 * t181 + t3 * t175 + t2 * t176 + t94 * t29 + (-t7 * t73 - t76 * t8) * mrSges(6,3) + (t100 * mrSges(6,3) + t102 * t180 + t104 * t185 + t106 * t183 + t59 * t109 + t11 * t176 + t76 * t188) * qJD(5) + t91) * t74 + (t83 + (t81 + t90 + t198) * t77) * qJD(4) + (((Ifges(5,5) * t74 + Ifges(5,6) * t77) * t202 + t118 + Ifges(4,4) * t199 + ((mrSges(5,2) * qJD(4) + t43) * t77 + (mrSges(5,1) * qJD(4) + t143) * t74 - t89) * qJ(2) - t86) * t75 + (t119 - t69 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t138 + t87 * t77 + t82 * t74 - t191) * t78) * qJD(1); (t59 * t108 + t81) * qJD(5) + (-mrSges(5,1) - t109) * t29 + t110 * mrSges(6,3) + (-m(6) * (t100 + t60) + t93) * t59 + t66 + t104 * t186 + t106 * t187 + t73 * t189 + t2 * t175 + t102 * t181 + (-t143 + t172) * t60 + (-t196 - t198 - t166 / 0.2e1 + t103 * t180 + t107 * t183 + t105 * t185 - t190) * t61 + (Ifges(6,5) * t183 + Ifges(6,6) * t185 + Ifges(6,3) * t180 + t193 - t197 + t203) * t62 - (Ifges(5,1) * t61 - t170 + t9) * t62 / 0.2e1 - (-Ifges(5,2) * t62 + t147 + t27 + t55) * t61 / 0.2e1 - t195; -t59 * (mrSges(6,1) * t35 + mrSges(6,2) * t34) + (Ifges(6,1) * t34 - t169) * t183 + t10 * t182 + (Ifges(6,5) * t34 - Ifges(6,6) * t35) * t180 - t36 * t23 + t37 * t24 + (t34 * t36 + t35 * t37) * mrSges(6,3) + t111 + t1 + (-Ifges(6,2) * t35 + t11 + t33) * t185;];
tauc = t4(:);
