% Calculate time derivative of joint inertia matrix for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:10
% EndTime: 2019-03-08 23:11:16
% DurationCPUTime: 2.84s
% Computational Cost: add. (3119->321), mult. (7583->457), div. (0->0), fcn. (7001->10), ass. (0->156)
t239 = mrSges(5,3) + mrSges(6,1);
t238 = mrSges(6,2) - mrSges(5,1);
t107 = sin(qJ(6));
t109 = cos(qJ(6));
t178 = qJD(6) * t109;
t230 = qJD(3) + qJD(4);
t108 = sin(qJ(4));
t210 = sin(qJ(3));
t212 = cos(qJ(4));
t213 = cos(qJ(3));
t84 = t108 * t213 + t212 * t210;
t64 = t230 * t84;
t150 = t212 * t213;
t83 = t108 * t210 - t150;
t137 = t107 * t64 + t83 * t178;
t179 = qJD(6) * t107;
t136 = -t109 * t64 + t83 * t179;
t196 = Ifges(7,4) * t109;
t91 = -Ifges(7,2) * t107 + t196;
t197 = Ifges(7,4) * t107;
t92 = Ifges(7,1) * t109 - t197;
t237 = t107 * t92 + t109 * t91;
t106 = sin(pkin(6));
t214 = cos(qJ(2));
t161 = t106 * t214;
t211 = sin(qJ(2));
t160 = t106 * t211;
t184 = cos(pkin(6));
t124 = t210 * t160 - t184 * t213;
t118 = t212 * t124;
t77 = t213 * t160 + t184 * t210;
t47 = t108 * t77 + t118;
t131 = -t107 * t47 + t109 * t161;
t236 = qJD(6) * t131;
t181 = t107 ^ 2 + t109 ^ 2;
t235 = (-mrSges(5,2) * t212 + (-t181 * mrSges(7,3) + t238) * t108) * pkin(3) * qJD(4);
t234 = -pkin(9) - pkin(8);
t198 = mrSges(5,2) - mrSges(6,3);
t90 = mrSges(7,1) * t107 + mrSges(7,2) * t109;
t233 = mrSges(6,3) + t90;
t232 = Ifges(4,1) - Ifges(4,2);
t101 = -t213 * pkin(3) - pkin(2);
t133 = -t84 * qJ(5) + t101;
t215 = pkin(4) + pkin(10);
t41 = t215 * t83 + t133;
t132 = t234 * t210;
t85 = t212 * t132;
t93 = t234 * t213;
t72 = -t108 * t93 - t85;
t45 = pkin(5) * t84 + t72;
t19 = -t107 * t41 + t109 * t45;
t156 = t210 * qJD(3);
t102 = pkin(3) * t156;
t63 = (t210 * qJD(4) + t156) * t108 - t230 * t150;
t135 = qJ(5) * t63 - qJD(5) * t84 + t102;
t11 = t215 * t64 + t135;
t20 = t107 * t45 + t109 * t41;
t183 = qJD(6) * t20;
t128 = qJD(3) * t93;
t129 = t108 * t132;
t157 = qJD(4) * t212;
t31 = -t212 * t128 + t230 * t129 - t93 * t157;
t22 = -t63 * pkin(5) + t31;
t2 = -t107 * t11 + t109 * t22 - t183;
t1 = qJD(6) * t19 + t107 * t22 + t109 * t11;
t208 = t1 * t107;
t23 = mrSges(7,2) * t63 - mrSges(7,3) * t136;
t24 = -mrSges(7,1) * t63 - mrSges(7,3) * t137;
t115 = m(7) * (t208 + t2 * t109 + (-t107 * t19 + t109 * t20) * qJD(6)) + t109 * t24 + t107 * t23;
t191 = t107 * t83;
t56 = mrSges(7,1) * t84 - mrSges(7,3) * t191;
t186 = t109 * t83;
t57 = -mrSges(7,2) * t84 + mrSges(7,3) * t186;
t229 = t57 * t178 - t56 * t179 + t115;
t228 = t84 * t239;
t205 = t63 * mrSges(6,1);
t227 = m(6) * t31 - t205;
t16 = mrSges(7,1) * t136 + mrSges(7,2) * t137;
t180 = qJD(4) * t108;
t30 = -t108 * t128 - t93 * t180 - t230 * t85;
t21 = -pkin(5) * t64 - t30;
t224 = m(7) * t21 - t64 * mrSges(6,1) + t16;
t73 = -t212 * t93 + t129;
t46 = -t83 * pkin(5) + t73;
t144 = mrSges(7,1) * t109 - mrSges(7,2) * t107;
t49 = t144 * t83;
t223 = -m(7) * t46 + t83 * mrSges(6,1) + t49;
t217 = m(5) * pkin(3);
t221 = -t108 * t217 + t198;
t220 = 0.2e1 * m(6);
t219 = 0.2e1 * m(7);
t86 = t144 * qJD(6);
t218 = 0.2e1 * t86;
t209 = pkin(3) * t108;
t159 = qJD(2) * t211;
t148 = t106 * t159;
t149 = qJD(2) * t161;
t113 = t77 * qJD(3) + t210 * t149;
t120 = t124 * qJD(3);
t114 = t213 * t149 - t120;
t121 = t108 * t124;
t15 = -qJD(4) * t121 + t108 * t114 + t212 * t113 + t77 * t157;
t36 = t107 * t161 + t109 * t47;
t4 = qJD(6) * t36 + t107 * t15 + t109 * t148;
t207 = t107 * t4;
t14 = qJD(4) * t118 + t108 * t113 - t212 * t114 + t77 * t180;
t48 = t212 * t77 - t121;
t5 = t48 * t14;
t204 = t63 * mrSges(5,3);
t202 = t64 * mrSges(5,3);
t200 = t83 * mrSges(5,3);
t152 = pkin(3) * t157;
t94 = t152 + qJD(5);
t98 = qJ(5) + t209;
t199 = t94 * t98;
t194 = t107 * mrSges(7,3);
t189 = t109 * mrSges(7,3);
t173 = t212 * pkin(3);
t172 = pkin(3) * t180;
t168 = mrSges(7,3) * t179;
t165 = mrSges(7,3) * t178;
t162 = t186 / 0.2e1;
t158 = qJD(3) * t213;
t100 = -t173 - pkin(4);
t147 = -t14 * t90 + t238 * t15 + t36 * t168 + t48 * t86;
t146 = -t14 * t98 + t48 * t94;
t145 = -t73 * t30 + t72 * t31;
t143 = Ifges(7,1) * t107 + t196;
t142 = Ifges(7,2) * t109 + t197;
t141 = Ifges(7,5) * t107 + Ifges(7,6) * t109;
t140 = -qJ(5) * t14 + t48 * qJD(5);
t139 = qJ(5) * t94 + qJD(5) * t98;
t138 = t181 * t172;
t134 = t106 ^ 2 * t214 * t159;
t130 = t137 * Ifges(7,5) - t136 * Ifges(7,6) - Ifges(7,3) * t63;
t127 = -t73 * t14 + t72 * t15 - t30 * t48 + t31 * t47;
t88 = t142 * qJD(6);
t89 = t143 * qJD(6);
t125 = -t237 * qJD(6) + t107 * t88 - t109 * t89;
t3 = -t107 * t148 + t109 * t15 + t236;
t122 = t207 + t109 * t3 + (-t107 * t36 - t109 * t131) * qJD(6);
t119 = m(7) * t122;
t39 = Ifges(7,6) * t84 + t142 * t83;
t40 = Ifges(7,5) * t84 + t143 * t83;
t8 = Ifges(7,4) * t137 - Ifges(7,2) * t136 - t63 * Ifges(7,6);
t9 = Ifges(7,1) * t137 - Ifges(7,4) * t136 - t63 * Ifges(7,5);
t116 = t19 * t168 + t21 * t90 - t39 * t178 / 0.2e1 + t46 * t86 - t107 * t8 / 0.2e1 - t89 * t191 / 0.2e1 - t88 * t162 + t109 * t9 / 0.2e1 + (-Ifges(7,5) * t109 / 0.2e1 + Ifges(7,6) * t107 / 0.2e1 + Ifges(6,4) - Ifges(5,5)) * t63 + t238 * t31 - (t83 * t91 + t40) * t179 / 0.2e1 + (t92 * t162 - t84 * t141 / 0.2e1) * qJD(6) + (Ifges(6,5) - Ifges(5,6) + t237 / 0.2e1) * t64;
t97 = -pkin(10) + t100;
t87 = (mrSges(4,1) * t210 + mrSges(4,2) * t213) * qJD(3);
t68 = -mrSges(6,2) * t83 - mrSges(6,3) * t84;
t67 = mrSges(5,1) * t83 + mrSges(5,2) * t84;
t50 = t83 * pkin(4) + t133;
t27 = mrSges(5,1) * t64 - mrSges(5,2) * t63;
t26 = -mrSges(6,2) * t64 + mrSges(6,3) * t63;
t25 = pkin(4) * t64 + t135;
t6 = [0.2e1 * m(7) * (-t131 * t4 + t3 * t36 - t5) + 0.2e1 * m(4) * (t113 * t124 + t114 * t77 - t134) + 0.2e1 * (m(6) + m(5)) * (t47 * t15 - t134 - t5); m(5) * ((t101 * t159 - t214 * t102) * t106 + t127) + m(6) * ((t50 * t159 - t214 * t25) * t106 + t127) - mrSges(3,2) * t149 + m(7) * (-t1 * t131 + t19 * t3 + t2 * t36 + t20 * t4) + t36 * t24 - t131 * t23 + t3 * t56 + t4 * t57 + (-t204 - t205) * t47 + t15 * t228 + (-t202 + t224) * t48 + (-t87 - t27 - t26) * t161 + (t200 + t223) * t14 + (-m(4) * pkin(2) - t213 * mrSges(4,1) + t210 * mrSges(4,2) - mrSges(3,1) + t67 + t68) * t148 + (m(4) * pkin(8) + mrSges(4,3)) * (t113 * t210 - t77 * t156 + (t114 + t120) * t213); 0.2e1 * t67 * t102 + 0.2e1 * m(5) * (t101 * t102 + t145) + t84 * t130 + 0.2e1 * ((-Ifges(5,4) - Ifges(6,6)) * t84 + (Ifges(6,3) + Ifges(5,2)) * t83) * t64 + (t1 * t20 + t19 * t2 + t21 * t46) * t219 + (t50 * t25 + t145) * t220 + (-0.2e1 * Ifges(4,4) * t210 + t232 * t213) * t156 + (0.2e1 * Ifges(4,4) * t213 + t232 * t210) * t158 + t137 * t40 - t136 * t39 + ((-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t84 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) - t141) * t83) * t63 + t8 * t186 + t9 * t191 + 0.2e1 * t20 * t23 + 0.2e1 * t19 * t24 + 0.2e1 * t46 * t16 - 0.2e1 * t21 * t49 + 0.2e1 * t50 * t26 + 0.2e1 * t2 * t56 + 0.2e1 * t1 * t57 + 0.2e1 * t25 * t68 - 0.2e1 * pkin(2) * t87 + 0.2e1 * t101 * t27 + 0.2e1 * t239 * (t30 * t83 + t31 * t84 - t63 * t72 - t64 * t73); m(6) * (t100 * t15 + t47 * t172 + t146) + m(7) * ((-t107 * t131 + t109 * t36) * t172 + t122 * t97 + t146) - t4 * t194 + t131 * t165 - t3 * t189 + (-t212 * t15 + (t108 * t47 + t212 * t48) * qJD(4)) * t217 - t113 * mrSges(4,1) - t114 * mrSges(4,2) + t147 + t221 * t14; -t20 * t165 - t202 * t209 - t152 * t200 - t1 * t194 - t2 * t189 - Ifges(4,6) * t156 + (-t212 * t31 + (t108 * t72 + t212 * t73) * qJD(4)) * t217 + t116 + Ifges(4,5) * t158 + t173 * t204 + t227 * t100 + t224 * t98 + (m(6) * t73 - t223) * t94 + (-m(6) * t98 + t221) * t30 + (-mrSges(4,1) * t158 + mrSges(4,2) * t156) * pkin(8) + t229 * t97 + (t109 * t56 + t107 * t57 + m(7) * (t107 * t20 + t109 * t19) + m(6) * t72 + t228) * t172; t98 * t218 + 0.2e1 * t233 * t94 + 0.2e1 * t235 + (t138 * t97 + t199) * t219 + (t100 * t172 + t199) * t220 + t125; t198 * t14 + m(6) * (-pkin(4) * t15 + t140) + m(7) * t140 - t215 * t119 + (-t207 + (-t3 + t236) * t109) * mrSges(7,3) + t147; -t229 * t215 + (-t208 + (-t2 - t183) * t109) * mrSges(7,3) + t198 * t30 + m(7) * (qJ(5) * t21 + qJD(5) * t46) + m(6) * (-pkin(4) * t31 - qJ(5) * t30 + qJD(5) * t73) + (pkin(4) * t63 - qJ(5) * t64 - qJD(5) * t83) * mrSges(6,1) + t116 + qJ(5) * t16 - qJD(5) * t49; (qJ(5) + t98) * t86 + t235 + m(7) * (-t138 * t215 + t139) + m(6) * (-pkin(4) * t172 + t139) + t125 + t233 * (qJD(5) + t94); qJ(5) * t218 + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t233) * qJD(5) + t125; m(6) * t15 + t119; (-t107 * t56 + t109 * t57) * qJD(6) + t115 + t227; (m(7) * t181 + m(6)) * t172; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t130; t144 * t172 + ((-mrSges(7,2) * t97 - Ifges(7,6)) * t109 + (-mrSges(7,1) * t97 - Ifges(7,5)) * t107) * qJD(6); ((mrSges(7,2) * t215 - Ifges(7,6)) * t109 + (mrSges(7,1) * t215 - Ifges(7,5)) * t107) * qJD(6); -t90 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
