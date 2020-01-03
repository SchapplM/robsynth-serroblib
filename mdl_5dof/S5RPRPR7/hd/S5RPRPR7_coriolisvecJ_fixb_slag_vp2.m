% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:54
% EndTime: 2019-12-31 18:19:01
% DurationCPUTime: 2.87s
% Computational Cost: add. (2525->325), mult. (6282->465), div. (0->0), fcn. (4042->8), ass. (0->166)
t102 = cos(qJ(3));
t154 = cos(pkin(9));
t100 = sin(qJ(3));
t96 = sin(pkin(9));
t156 = t96 * t100;
t109 = t102 * t154 - t156;
t73 = t109 * qJD(1);
t202 = -t73 / 0.2e1;
t201 = -qJD(3) / 0.2e1;
t101 = cos(qJ(5));
t91 = sin(pkin(8)) * pkin(1) + pkin(6);
t83 = t91 * qJD(1);
t124 = qJ(4) * qJD(1) + t83;
t141 = t100 * qJD(2);
t57 = t102 * t124 + t141;
t49 = t154 * t57;
t95 = t102 * qJD(2);
t56 = -t100 * t124 + t95;
t51 = qJD(3) * pkin(3) + t56;
t24 = t96 * t51 + t49;
t22 = qJD(3) * pkin(7) + t24;
t136 = -cos(pkin(8)) * pkin(1) - pkin(2);
t82 = -pkin(3) * t102 + t136;
t149 = qJD(1) * t82;
t72 = qJD(4) + t149;
t126 = t154 * t100;
t145 = qJD(1) * t102;
t74 = -qJD(1) * t126 - t145 * t96;
t30 = -t73 * pkin(4) + t74 * pkin(7) + t72;
t99 = sin(qJ(5));
t7 = t101 * t30 - t22 * t99;
t200 = t7 * mrSges(6,1);
t169 = t74 * mrSges(5,3);
t54 = qJD(3) * t101 + t74 * t99;
t55 = qJD(3) * t99 - t101 * t74;
t162 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t54 + mrSges(6,2) * t55 - t169;
t139 = qJD(1) * qJD(3);
t129 = t100 * t139;
t199 = qJD(5) - t73;
t140 = t100 * qJD(4);
t198 = -qJD(1) * t140 - t57 * qJD(3);
t76 = t109 * qJD(3);
t66 = qJD(1) * t76;
t33 = qJD(5) * t54 + t101 * t66;
t81 = t96 * t102 + t126;
t75 = t81 * qJD(3);
t65 = qJD(1) * t75;
t19 = mrSges(6,1) * t65 - mrSges(6,3) * t33;
t34 = -qJD(5) * t55 - t66 * t99;
t20 = -mrSges(6,2) * t65 + mrSges(6,3) * t34;
t197 = t101 * t20 - t99 * t19;
t146 = qJD(1) * t100;
t143 = qJD(4) * t102;
t144 = qJD(3) * t100;
t61 = qJD(3) * t95 - t144 * t83;
t47 = (-qJ(4) * t144 + t143) * qJD(1) + t61;
t15 = t154 * t47 + t198 * t96;
t123 = pkin(3) * t129;
t35 = pkin(4) * t65 - pkin(7) * t66 + t123;
t1 = qJD(5) * t7 + t101 * t15 + t35 * t99;
t8 = t101 * t22 + t30 * t99;
t2 = -qJD(5) * t8 + t101 * t35 - t15 * t99;
t196 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t33 + Ifges(6,6) * t34;
t167 = t74 * Ifges(5,4);
t192 = Ifges(5,6) * t201 + t167 / 0.2e1 + Ifges(5,2) * t202;
t118 = mrSges(6,1) * t99 + mrSges(6,2) * t101;
t166 = t96 * t57;
t23 = t154 * t51 - t166;
t21 = -qJD(3) * pkin(4) - t23;
t110 = t21 * t118;
t111 = Ifges(6,5) * t101 - Ifges(6,6) * t99;
t160 = Ifges(6,4) * t101;
t113 = -Ifges(6,2) * t99 + t160;
t180 = Ifges(6,4) * t99;
t115 = Ifges(6,1) * t101 - t180;
t181 = Ifges(6,4) * t55;
t17 = Ifges(6,2) * t54 + Ifges(6,6) * t199 + t181;
t164 = t99 * t17;
t52 = Ifges(6,4) * t54;
t18 = Ifges(6,1) * t55 + Ifges(6,5) * t199 + t52;
t182 = t101 / 0.2e1;
t189 = t55 / 0.2e1;
t195 = -t164 / 0.2e1 + t18 * t182 + t199 * t111 / 0.2e1 + t115 * t189 + t54 * t113 / 0.2e1 + t110;
t194 = t33 / 0.2e1;
t193 = t34 / 0.2e1;
t191 = -t54 / 0.2e1;
t190 = -t55 / 0.2e1;
t188 = t65 / 0.2e1;
t187 = -t199 / 0.2e1;
t184 = -t99 / 0.2e1;
t183 = pkin(3) * t96;
t70 = Ifges(5,4) * t73;
t179 = t1 * t101;
t14 = -t154 * t198 + t47 * t96;
t155 = qJ(4) + t91;
t79 = t155 * t102;
t44 = t126 * t155 + t79 * t96;
t178 = t14 * t44;
t177 = t14 * t109;
t176 = t54 * Ifges(6,6);
t175 = t55 * Ifges(6,5);
t174 = t65 * mrSges(5,3);
t173 = t66 * mrSges(5,3);
t172 = t199 * Ifges(6,3);
t171 = t73 * mrSges(5,3);
t168 = t74 * Ifges(5,1);
t165 = t99 * mrSges(6,3);
t161 = Ifges(4,4) * t100;
t157 = t101 * t73;
t153 = Ifges(4,5) * qJD(3);
t152 = Ifges(5,5) * qJD(3);
t151 = Ifges(4,6) * qJD(3);
t148 = qJD(5) * t99;
t142 = qJD(5) * t101;
t138 = pkin(3) * t146;
t135 = mrSges(4,3) * t146;
t134 = mrSges(4,3) * t145;
t94 = Ifges(4,4) * t145;
t132 = m(4) * t91 + mrSges(4,3);
t131 = t154 * pkin(3);
t130 = t65 * mrSges(5,1) + t66 * mrSges(5,2);
t125 = qJD(3) * t155;
t121 = -t1 * t99 - t101 * t2;
t120 = t101 * t8 - t7 * t99;
t119 = -t7 * t101 - t8 * t99;
t85 = t136 * qJD(1);
t117 = mrSges(6,1) * t101 - mrSges(6,2) * t99;
t116 = Ifges(6,1) * t99 + t160;
t114 = Ifges(6,2) * t101 + t180;
t112 = Ifges(6,5) * t99 + Ifges(6,6) * t101;
t40 = -pkin(4) * t109 - t81 * pkin(7) + t82;
t45 = t154 * t79 - t155 * t156;
t13 = t101 * t45 + t40 * t99;
t12 = t101 * t40 - t45 * t99;
t68 = t102 * t83 + t141;
t105 = -t102 * t125 - t140;
t103 = qJD(5) * t119 - t2 * t99 + t179;
t92 = -t131 - pkin(4);
t86 = -qJD(3) * mrSges(4,2) + t134;
t84 = qJD(3) * mrSges(4,1) - t135;
t78 = Ifges(4,1) * t146 + t153 + t94;
t77 = t151 + (t102 * Ifges(4,2) + t161) * qJD(1);
t67 = -t100 * t83 + t95;
t64 = Ifges(6,3) * t65;
t62 = t68 * qJD(3);
t59 = -qJD(3) * mrSges(5,2) + t171;
t58 = -t100 * t125 + t143;
t46 = -mrSges(5,1) * t73 - mrSges(5,2) * t74;
t42 = t152 + t70 - t168;
t39 = pkin(3) * t144 + pkin(4) * t75 - pkin(7) * t76;
t38 = -pkin(4) * t74 - pkin(7) * t73 + t138;
t37 = mrSges(6,1) * t199 - mrSges(6,3) * t55;
t36 = -mrSges(6,2) * t199 + mrSges(6,3) * t54;
t29 = t105 * t96 + t154 * t58;
t28 = -t105 * t154 + t58 * t96;
t26 = t154 * t56 - t166;
t25 = t56 * t96 + t49;
t16 = t172 + t175 + t176;
t11 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t10 = t101 * t26 + t38 * t99;
t9 = t101 * t38 - t26 * t99;
t6 = t33 * Ifges(6,1) + t34 * Ifges(6,4) + t65 * Ifges(6,5);
t5 = t33 * Ifges(6,4) + t34 * Ifges(6,2) + t65 * Ifges(6,6);
t4 = -qJD(5) * t13 + t101 * t39 - t29 * t99;
t3 = qJD(5) * t12 + t101 * t29 + t39 * t99;
t27 = [t82 * t130 + t29 * t59 + t44 * t11 + t3 * t36 + t4 * t37 + t12 * t19 + t13 * t20 + (t172 / 0.2e1 + t176 / 0.2e1 + t175 / 0.2e1 - t8 * mrSges(6,2) + t200 + t72 * mrSges(5,1) + t16 / 0.2e1 + 0.2e1 * t192) * t75 + (-t168 / 0.2e1 + t72 * mrSges(5,2) + t70 / 0.2e1 + t152 / 0.2e1 + t42 / 0.2e1 + t119 * mrSges(6,3) + t195) * t76 + t162 * t28 + m(6) * (t1 * t13 + t12 * t2 + t21 * t28 + t3 * t8 + t4 * t7 + t178) + m(5) * (t15 * t45 - t23 * t28 + t24 * t29 + t178) + (-t23 * t76 - t24 * t75 + t44 * t66 - t45 * t65) * mrSges(5,3) + (t132 * t61 + (t78 / 0.2e1 - t91 * t84 + 0.3e1 / 0.2e1 * t94 + t153 / 0.2e1 - t132 * t67 + 0.2e1 * t85 * mrSges(4,2)) * qJD(3)) * t102 - (t64 / 0.2e1 - Ifges(5,4) * t66 - t15 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t65 + t196) * t109 + (t111 * t188 + t115 * t194 + t113 * t193 - Ifges(5,4) * t65 + Ifges(5,1) * t66 + t6 * t182 + t5 * t184 + (mrSges(5,3) + t118) * t14 + t121 * mrSges(6,3) + (t21 * t117 + t112 * t187 + t114 * t191 + t116 * t190 + t18 * t184 - t101 * t17 / 0.2e1 - t120 * mrSges(6,3)) * qJD(5)) * t81 + (t132 * t62 + (t85 * mrSges(4,1) - t77 / 0.2e1 - t91 * t86 - t151 / 0.2e1 - t132 * t68 + (t136 * mrSges(4,1) - 0.3e1 / 0.2e1 * t161 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t102) * qJD(1) + (m(5) * (t72 + t149) + t46 + qJD(1) * (-mrSges(5,1) * t109 + mrSges(5,2) * t81)) * pkin(3)) * qJD(3)) * t100; -(t11 + t173) * t109 + t162 * t75 + (t101 * t36 - t99 * t37 + t59) * t76 + (-t100 * t84 + t102 * t86 + (-t100 ^ 2 - t102 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t174 + (-t101 * t37 - t99 * t36) * qJD(5) + t197) * t81 + m(4) * (t61 * t100 - t102 * t62 + (-t100 * t67 + t102 * t68) * qJD(3)) + m(5) * (t15 * t81 - t23 * t75 + t24 * t76 - t177) + m(6) * (t103 * t81 + t120 * t76 + t21 * t75 - t177); (-t148 * t8 + t179 + (-t142 + t157) * t7) * mrSges(6,3) + t99 * t6 / 0.2e1 + t92 * t11 + (t135 + t84) * t68 - (-Ifges(4,2) * t146 + t78 + t94) * t145 / 0.2e1 + (Ifges(5,1) * t73 + t16 + t167) * t74 / 0.2e1 + ((-t14 * t154 + t15 * t96) * pkin(3) - t138 * t72 + t23 * t25 - t24 * t26) * m(5) + (-m(6) * t21 - t162) * t25 + (-t85 * (mrSges(4,1) * t100 + mrSges(4,2) * t102) - (Ifges(4,1) * t102 - t161) * t146 / 0.2e1) * qJD(1) + (m(6) * t103 - t142 * t37 - t148 * t36 + t197) * (pkin(7) + t183) - t72 * (-mrSges(5,1) * t74 + mrSges(5,2) * t73) - Ifges(5,6) * t65 + Ifges(5,5) * t66 - t26 * t59 - t61 * mrSges(4,2) - t62 * mrSges(4,1) - t10 * t36 - t9 * t37 - t15 * mrSges(5,2) - t18 * t157 / 0.2e1 + t77 * t146 / 0.2e1 - t46 * t138 + (m(6) * t92 - mrSges(5,1) - t117) * t14 - Ifges(4,6) * t129 / 0.2e1 - t174 * t183 - t131 * t173 + t195 * qJD(5) + t73 * t164 / 0.2e1 - t8 * (mrSges(6,2) * t74 - t165 * t73) - t2 * t165 - t24 * t169 + (t134 - t86) * t67 + (Ifges(5,5) * t73 + Ifges(5,6) * t74) * t201 + (Ifges(5,2) * t74 + t42 + t70) * t202 + t74 * t200 + t116 * t194 + (-Ifges(6,3) * t74 + t111 * t73) * t187 + t112 * t188 + (-Ifges(6,5) * t74 + t115 * t73) * t190 + (-Ifges(6,6) * t74 + t113 * t73) * t191 + t74 * t192 + t114 * t193 + t5 * t182 - t73 * t110 - m(6) * (t10 * t8 + t7 * t9) + t23 * t171 + t139 * Ifges(4,5) * t102 / 0.2e1; -t73 * t59 + t162 * t74 + (-t199 * t37 + t20) * t99 + (t199 * t36 + t19) * t101 + t130 + (t199 * t120 + t21 * t74 - t121) * m(6) + (-t23 * t74 - t24 * t73 + t123) * m(5); t64 - t21 * (mrSges(6,1) * t55 + mrSges(6,2) * t54) + (Ifges(6,1) * t54 - t181) * t190 + t17 * t189 + (Ifges(6,5) * t54 - Ifges(6,6) * t55) * t187 - t7 * t36 + t8 * t37 + (t54 * t7 + t55 * t8) * mrSges(6,3) + (-Ifges(6,2) * t55 + t18 + t52) * t191 + t196;];
tauc = t27(:);
