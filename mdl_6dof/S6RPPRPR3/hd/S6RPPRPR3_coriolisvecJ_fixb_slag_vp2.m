% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:40:19
% EndTime: 2018-11-23 15:40:22
% DurationCPUTime: 3.20s
% Computational Cost: add. (3499->370), mult. (7507->519), div. (0->0), fcn. (4624->8), ass. (0->174)
t111 = cos(qJ(6));
t197 = t111 / 0.2e1;
t216 = m(5) + m(4);
t109 = sin(qJ(6));
t176 = Ifges(7,4) * t109;
t215 = Ifges(7,2) * t197 + t176 / 0.2e1;
t106 = sin(pkin(10));
t112 = cos(qJ(4));
t169 = cos(pkin(10));
t138 = t169 * t112;
t134 = qJD(1) * t138;
t110 = sin(qJ(4));
t162 = qJD(1) * t110;
t75 = t106 * t162 - t134;
t185 = t75 * mrSges(6,3);
t58 = qJD(4) * t111 + t109 * t75;
t59 = qJD(4) * t109 - t111 * t75;
t180 = -qJD(4) * mrSges(6,1) - mrSges(7,1) * t58 + mrSges(7,2) * t59 - t185;
t127 = mrSges(7,1) * t109 + mrSges(7,2) * t111;
t153 = t112 * qJD(2);
t95 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t82 = qJD(1) * t95 + qJD(3);
t209 = t110 * (qJ(5) * qJD(1) - t82) - t153;
t174 = t106 * t209;
t161 = qJD(1) * t112;
t155 = t110 * qJD(2);
t61 = t112 * t82 - t155;
t52 = -qJ(5) * t161 + t61;
t51 = qJD(4) * pkin(4) + t52;
t23 = t169 * t51 + t174;
t21 = -qJD(4) * pkin(5) - t23;
t214 = t21 * t127;
t97 = sin(pkin(9)) * pkin(1) + qJ(3);
t88 = qJD(1) * t97;
t117 = -t106 * t112 - t110 * t169;
t76 = t117 * qJD(1);
t213 = qJD(6) - t76;
t158 = qJD(5) * t112;
t212 = -qJD(1) * t158 + t209 * qJD(4);
t69 = qJD(4) * t76;
t34 = qJD(6) * t58 + t111 * t69;
t152 = qJD(1) * qJD(4);
t142 = t110 * t152;
t68 = -qJD(4) * t134 + t106 * t142;
t19 = -mrSges(7,1) * t68 - mrSges(7,3) * t34;
t35 = -qJD(6) * t59 - t109 * t69;
t20 = mrSges(7,2) * t68 + mrSges(7,3) * t35;
t211 = t109 * t19 - t111 * t20;
t154 = t110 * qJD(5);
t159 = qJD(4) * t112;
t56 = -qJD(4) * t155 + t82 * t159;
t46 = (-qJ(5) * t159 - t154) * qJD(1) + t56;
t15 = t106 * t212 + t169 * t46;
t105 = qJD(1) * qJD(3);
t141 = t112 * t152;
t86 = pkin(4) * t141 + t105;
t31 = -pkin(5) * t68 - pkin(8) * t69 + t86;
t49 = t169 * t209;
t24 = t106 * t51 - t49;
t22 = qJD(4) * pkin(8) + t24;
t74 = pkin(4) * t162 + qJD(5) + t88;
t30 = -pkin(5) * t76 + pkin(8) * t75 + t74;
t5 = -t109 * t22 + t111 * t30;
t1 = qJD(6) * t5 + t109 * t31 + t111 * t15;
t6 = t109 * t30 + t111 * t22;
t2 = -qJD(6) * t6 - t109 * t15 + t111 * t31;
t210 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t34 + Ifges(7,6) * t35;
t36 = -mrSges(7,2) * t213 + mrSges(7,3) * t58;
t37 = mrSges(7,1) * t213 - mrSges(7,3) * t59;
t120 = t109 * t36 + t111 * t37;
t188 = t68 * mrSges(6,3);
t208 = qJD(6) * t120 - t188 + t211;
t131 = t6 * t109 + t5 * t111;
t207 = qJD(6) * t131 - t1 * t111 + t109 * t2;
t205 = -t58 / 0.2e1;
t204 = -t59 / 0.2e1;
t203 = t59 / 0.2e1;
t202 = -t213 / 0.2e1;
t199 = t109 / 0.2e1;
t198 = -t111 / 0.2e1;
t196 = Ifges(6,4) * t75;
t71 = Ifges(6,4) * t76;
t195 = Ifges(7,4) * t59;
t194 = pkin(4) * t106;
t14 = t106 * t46 - t169 * t212;
t170 = qJ(5) - t95;
t79 = t170 * t110;
t44 = -t106 * t79 + t138 * t170;
t193 = t14 * t44;
t83 = t106 * t110 - t138;
t192 = t14 * t83;
t191 = t14 * t117;
t190 = t58 * Ifges(7,6);
t189 = t59 * Ifges(7,5);
t187 = t69 * mrSges(6,3);
t186 = t213 * Ifges(7,3);
t184 = t75 * Ifges(6,1);
t183 = t76 * mrSges(6,3);
t182 = t76 * Ifges(6,2);
t179 = mrSges(7,3) * t111;
t178 = Ifges(5,4) * t110;
t177 = Ifges(5,4) * t112;
t175 = Ifges(7,4) * t111;
t173 = t109 * mrSges(7,3);
t168 = mrSges(4,3) * qJD(1);
t167 = Ifges(5,5) * qJD(4);
t166 = Ifges(6,5) * qJD(4);
t165 = Ifges(5,6) * qJD(4);
t164 = Ifges(6,6) * qJD(4);
t160 = qJD(4) * t110;
t157 = qJD(6) * t109;
t156 = qJD(6) * t111;
t92 = pkin(4) * t159 + qJD(3);
t151 = pkin(4) * t161;
t150 = mrSges(5,3) * t162;
t149 = mrSges(5,3) * t161;
t17 = Ifges(7,2) * t58 + Ifges(7,6) * t213 + t195;
t148 = t17 * t199;
t55 = Ifges(7,4) * t58;
t18 = Ifges(7,1) * t59 + Ifges(7,5) * t213 + t55;
t147 = t18 * t198;
t87 = t110 * pkin(4) + t97;
t146 = m(5) * t95 - mrSges(5,3);
t145 = t169 * pkin(4);
t144 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t11 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t143 = t11 + t187;
t139 = t170 * t112;
t137 = 0.2e1 * t88;
t132 = t1 * t109 + t2 * t111;
t130 = t5 * t109 - t6 * t111;
t77 = -qJD(4) * t138 + t106 * t160;
t78 = t117 * qJD(4);
t129 = t23 * t78 - t24 * t77;
t128 = -mrSges(7,1) * t111 + mrSges(7,2) * t109;
t126 = Ifges(7,1) * t111 - t176;
t125 = Ifges(7,1) * t109 + t175;
t124 = -Ifges(7,2) * t109 + t175;
t122 = Ifges(7,5) * t111 - Ifges(7,6) * t109;
t121 = Ifges(7,5) * t109 + Ifges(7,6) * t111;
t41 = -pkin(5) * t117 + pkin(8) * t83 + t87;
t45 = -t106 * t139 - t169 * t79;
t12 = -t109 * t45 + t111 * t41;
t13 = t109 * t41 + t111 * t45;
t62 = t110 * t82 + t153;
t63 = -qJD(4) * mrSges(6,2) + t183;
t119 = t109 * t37 - t111 * t36 - t63;
t115 = t160 * t170 - t158;
t100 = -t145 - pkin(5);
t90 = qJD(4) * mrSges(5,1) - t149;
t89 = -qJD(4) * mrSges(5,2) - t150;
t85 = qJD(1) * (t110 * mrSges(5,1) + mrSges(5,2) * t112);
t81 = t167 + (t112 * Ifges(5,1) - t178) * qJD(1);
t80 = t165 + (-Ifges(5,2) * t110 + t177) * qJD(1);
t66 = Ifges(7,3) * t68;
t60 = -qJD(4) * t139 - t154;
t57 = t62 * qJD(4);
t47 = -mrSges(6,1) * t76 - mrSges(6,2) * t75;
t43 = t166 + t71 - t184;
t42 = t164 + t182 - t196;
t40 = -t75 * pkin(5) - t76 * pkin(8) + t151;
t38 = -pkin(5) * t77 - pkin(8) * t78 + t92;
t28 = t106 * t115 + t169 * t60;
t27 = t106 * t60 - t115 * t169;
t26 = t169 * t52 + t174;
t25 = t106 * t52 - t49;
t16 = t186 + t189 + t190;
t10 = t109 * t40 + t111 * t26;
t9 = -t109 * t26 + t111 * t40;
t8 = t34 * Ifges(7,1) + t35 * Ifges(7,4) - t68 * Ifges(7,5);
t7 = t34 * Ifges(7,4) + t35 * Ifges(7,2) - t68 * Ifges(7,6);
t4 = -qJD(6) * t13 - t109 * t28 + t111 * t38;
t3 = qJD(6) * t12 + t109 * t38 + t111 * t28;
t29 = [t92 * t47 + t87 * t144 + t28 * t63 + t44 * t11 + t3 * t36 + t4 * t37 + t13 * t20 + t12 * t19 + (t6 * mrSges(7,2) - t186 / 0.2e1 - t189 / 0.2e1 - t190 / 0.2e1 - t5 * mrSges(7,1) - t196 / 0.2e1 + t182 / 0.2e1 - t74 * mrSges(6,1) + t164 / 0.2e1 + t42 / 0.2e1 - t16 / 0.2e1) * t77 - (-t214 + t122 * t202 + t126 * t204 + t124 * t205 + t147 + t148 + t184 / 0.2e1 - t71 / 0.2e1 - t74 * mrSges(6,2) - t166 / 0.2e1 - t43 / 0.2e1 + t131 * mrSges(7,3)) * t78 + t180 * t27 + m(6) * (t15 * t45 - t23 * t27 + t24 * t28 + t74 * t92 + t86 * t87 + t193) + m(7) * (t1 * t13 + t12 * t2 + t21 * t27 + t3 * t6 + t4 * t5 + t193) + (t44 * t69 + t45 * t68 - t129) * mrSges(6,3) + (t137 * t216 + 0.2e1 * t168 + t85) * qJD(3) + (mrSges(5,2) * t105 - t146 * t57 + (t95 * t89 - t80 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4) * t161 - t165 / 0.2e1 + t146 * t62 + t137 * mrSges(5,1)) * qJD(4)) * t112 - (-t15 * mrSges(6,3) - t66 / 0.2e1 - Ifges(6,4) * t69 + t86 * mrSges(6,1) + (-Ifges(7,3) / 0.2e1 - Ifges(6,2)) * t68 + t210) * t117 + (t8 * t198 + t7 * t199 - t34 * t126 / 0.2e1 - t35 * t124 / 0.2e1 - Ifges(6,1) * t69 - t86 * mrSges(6,2) + (-mrSges(6,3) - t127) * t14 + t132 * mrSges(7,3) + (t17 * t197 + t21 * t128 + t213 * t121 / 0.2e1 + t125 * t203 + t58 * t215 + t18 * t199 - t130 * mrSges(7,3)) * qJD(6) + (t122 / 0.2e1 - Ifges(6,4)) * t68) * t83 + (mrSges(5,1) * t105 + t146 * t56 + (-t81 / 0.2e1 - t88 * mrSges(5,2) - t95 * t90 - t167 / 0.2e1 - t146 * t61 + (-t97 * mrSges(5,2) + 0.3e1 / 0.2e1 * t178 + (0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,1)) * t112) * qJD(1)) * qJD(4)) * t110; -t143 * t117 - t180 * t77 - t119 * t78 + (-t110 * t89 - t112 * t90 + (-t110 ^ 2 - t112 ^ 2) * qJD(1) * mrSges(5,3)) * qJD(4) + t208 * t83 + m(5) * (t57 * t110 + t112 * t56 + (-t110 * t62 - t112 * t61) * qJD(4)) + m(6) * (-t15 * t83 + t23 * t77 + t24 * t78 - t191) + m(7) * (-t130 * t78 + t207 * t83 - t21 * t77 - t191); t143 * t83 - t180 * t78 + (-t110 * t90 + t112 * t89) * qJD(4) + t119 * t77 + t208 * t117 + m(5) * (t56 * t110 - t112 * t57 + (-t110 * t61 + t112 * t62) * qJD(4)) + m(7) * (t117 * t207 + t130 * t77 - t21 * t78 + t192) + m(6) * (-t117 * t15 + t129 + t192) + (-m(6) * t74 - m(7) * t131 - t216 * t88 - t120 - t168 - t47 - t85) * qJD(1); (t90 + t149) * t62 + (-t89 - t150) * t61 + t100 * t11 - t74 * (-mrSges(6,1) * t75 + mrSges(6,2) * t76) - qJD(4) * (Ifges(6,5) * t76 + Ifges(6,6) * t75) / 0.2e1 - t75 * t42 / 0.2e1 + Ifges(6,5) * t69 - t26 * t63 - t56 * mrSges(5,2) - t57 * mrSges(5,1) - t10 * t36 - t9 * t37 - t15 * mrSges(6,2) - (Ifges(6,2) * t75 + t43 + t71) * t76 / 0.2e1 + (Ifges(6,1) * t76 + t16 + t196) * t75 / 0.2e1 - t24 * t185 - t145 * t187 + (m(7) * t100 - mrSges(6,1) + t128) * t14 - t5 * (-mrSges(7,1) * t75 - t179 * t76) - t6 * (mrSges(7,2) * t75 - t173 * t76) - t2 * t173 + t213 * t214 + t35 * t215 + (-t156 * t5 - t157 * t6) * mrSges(7,3) + (-m(7) * t207 - t156 * t37 - t157 * t36 - t211) * (pkin(8) + t194) + t76 * t147 + t76 * t148 + (-t151 * t74 + t23 * t25 - t24 * t26 + (t106 * t15 - t14 * t169) * pkin(4)) * m(6) + (-m(7) * t21 - t180) * t25 + (-t88 * (mrSges(5,1) * t112 - mrSges(5,2) * t110) + (-t112 * (-Ifges(5,1) * t110 - t177) / 0.2e1 + t110 * (-Ifges(5,2) * t112 - t178) / 0.2e1) * qJD(1)) * qJD(1) + t80 * t161 / 0.2e1 + t81 * t162 / 0.2e1 + t18 * t156 / 0.2e1 - t17 * t157 / 0.2e1 - m(7) * (t10 * t6 + t5 * t9) - (-Ifges(5,5) * t110 - Ifges(5,6) * t112) * t152 / 0.2e1 - t47 * t151 - Ifges(5,5) * t142 - Ifges(5,6) * t141 + t34 * t125 / 0.2e1 + t1 * t179 + t23 * t183 + t188 * t194 + t7 * t197 + t8 * t199 + (-Ifges(7,3) * t75 + t122 * t76) * t202 + (-Ifges(7,5) * t75 + t126 * t76) * t204 + (-Ifges(7,6) * t75 + t124 * t76) * t205 + (Ifges(6,6) - t121 / 0.2e1) * t68 + (t122 * t213 + t124 * t58 + t126 * t59) * qJD(6) / 0.2e1; -t76 * t63 + t180 * t75 + (t213 * t36 + t19) * t111 + (-t213 * t37 + t20) * t109 + t144 + (-t130 * t213 + t21 * t75 + t132) * m(7) + (-t23 * t75 - t24 * t76 + t86) * m(6); -t66 - t21 * (mrSges(7,1) * t59 + mrSges(7,2) * t58) + (Ifges(7,1) * t58 - t195) * t204 + t17 * t203 + (Ifges(7,5) * t58 - Ifges(7,6) * t59) * t202 - t5 * t36 + t6 * t37 + (t5 * t58 + t59 * t6) * mrSges(7,3) + (-Ifges(7,2) * t59 + t18 + t55) * t205 + t210;];
tauc  = t29(:);
