% Calculate time derivative of joint inertia matrix for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:16
% EndTime: 2019-03-09 03:51:22
% DurationCPUTime: 3.07s
% Computational Cost: add. (7960->380), mult. (17764->578), div. (0->0), fcn. (18407->10), ass. (0->152)
t152 = sin(pkin(11));
t191 = pkin(8) + qJ(4);
t141 = t191 * t152;
t154 = cos(pkin(11));
t143 = t191 * t154;
t157 = sin(qJ(5));
t160 = cos(qJ(5));
t114 = -t157 * t141 + t160 * t143;
t153 = sin(pkin(10));
t192 = pkin(7) + qJ(2);
t142 = t192 * t153;
t155 = cos(pkin(10));
t144 = t192 * t155;
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t203 = -t161 * t142 - t144 * t158;
t176 = t154 * t160;
t163 = t152 * t157 - t176;
t127 = t163 * qJD(5);
t140 = t153 * t161 + t158 * t155;
t130 = t140 * qJD(3);
t139 = t152 * t160 + t154 * t157;
t128 = t139 * qJD(5);
t138 = t153 * t158 - t161 * t155;
t129 = t138 * qJD(3);
t61 = -t140 * t128 + t163 * t129;
t62 = t140 * t127 + t139 * t129;
t201 = Ifges(6,5) * t61 + Ifges(6,6) * t62 + Ifges(6,3) * t130;
t200 = 2 * m(5);
t199 = 2 * m(6);
t198 = 2 * m(7);
t150 = t154 ^ 2;
t196 = m(7) * pkin(5);
t193 = pkin(5) * t128;
t178 = t140 * t154;
t169 = -pkin(2) * t155 - pkin(1);
t104 = pkin(3) * t138 - qJ(4) * t140 + t169;
t115 = -t158 * t142 + t144 * t161;
t74 = t154 * t104 - t115 * t152;
t49 = pkin(4) * t138 - pkin(8) * t178 + t74;
t179 = t140 * t152;
t75 = t152 * t104 + t154 * t115;
t57 = -pkin(8) * t179 + t75;
t26 = t157 * t49 + t160 * t57;
t156 = sin(qJ(6));
t159 = cos(qJ(6));
t105 = -t139 * t156 - t159 * t163;
t69 = t105 * qJD(6) - t127 * t159 - t128 * t156;
t106 = t139 * t159 - t156 * t163;
t70 = -t106 * qJD(6) + t127 * t156 - t128 * t159;
t190 = Ifges(7,5) * t69 + Ifges(7,6) * t70;
t85 = pkin(3) * t130 + qJ(4) * t129 - qJD(4) * t140;
t89 = -t138 * qJD(2) + t203 * qJD(3);
t48 = t152 * t85 + t154 * t89;
t189 = Ifges(5,4) * t152;
t188 = Ifges(5,4) * t154;
t90 = t140 * qJD(2) + t115 * qJD(3);
t187 = t203 * t90;
t186 = t128 * mrSges(6,1);
t185 = t130 * Ifges(5,5);
t184 = t130 * Ifges(5,6);
t183 = t152 * Ifges(5,2);
t182 = t157 * t57;
t181 = t129 * t152;
t180 = t129 * t154;
t94 = -mrSges(5,1) * t181 - mrSges(5,2) * t180;
t175 = -Ifges(6,5) * t127 - Ifges(6,6) * t128;
t173 = qJD(5) * t160;
t172 = qJD(6) * t156;
t171 = qJD(6) * t159;
t95 = t139 * t140;
t96 = t163 * t140;
t63 = t156 * t96 - t159 * t95;
t18 = t63 * qJD(6) + t156 * t62 + t159 * t61;
t64 = -t156 * t95 - t159 * t96;
t19 = -t64 * qJD(6) - t156 * t61 + t159 * t62;
t170 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t130;
t147 = -pkin(4) * t154 - pkin(3);
t31 = -t62 * mrSges(6,1) + t61 * mrSges(6,2);
t8 = -t19 * mrSges(7,1) + t18 * mrSges(7,2);
t34 = -t70 * mrSges(7,1) + t69 * mrSges(7,2);
t47 = -t152 * t89 + t154 * t85;
t25 = t160 * t49 - t182;
t168 = t130 * mrSges(4,1) - t129 * mrSges(4,2);
t112 = -t160 * t141 - t143 * t157;
t121 = t127 * mrSges(6,2);
t166 = -t121 + t34;
t91 = pkin(4) * t179 - t203;
t165 = Ifges(5,5) * t154 - Ifges(5,6) * t152;
t20 = pkin(5) * t138 + pkin(9) * t96 + t25;
t21 = -pkin(9) * t95 + t26;
t9 = -t156 * t21 + t159 * t20;
t10 = t156 * t20 + t159 * t21;
t92 = -pkin(9) * t139 + t112;
t93 = -pkin(9) * t163 + t114;
t55 = -t156 * t93 + t159 * t92;
t56 = t156 * t92 + t159 * t93;
t87 = -t141 * t173 + qJD(4) * t176 + (-qJD(4) * t152 - qJD(5) * t143) * t157;
t77 = -pkin(9) * t128 + t87;
t88 = -t139 * qJD(4) - t114 * qJD(5);
t78 = pkin(9) * t127 + t88;
t23 = t55 * qJD(6) + t156 * t78 + t159 * t77;
t24 = -t56 * qJD(6) - t156 * t77 + t159 * t78;
t164 = t24 * mrSges(7,1) - t23 * mrSges(7,2) + t190;
t39 = pkin(4) * t130 + pkin(8) * t180 + t47;
t41 = pkin(8) * t181 + t48;
t12 = -t26 * qJD(5) - t157 * t41 + t160 * t39;
t4 = pkin(5) * t130 - pkin(9) * t61 + t12;
t11 = -qJD(5) * t182 + t157 * t39 + t160 * t41 + t49 * t173;
t5 = pkin(9) * t62 + t11;
t2 = t9 * qJD(6) + t156 * t4 + t159 * t5;
t3 = -t10 * qJD(6) - t156 * t5 + t159 * t4;
t162 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t170;
t76 = -pkin(4) * t181 + t90;
t132 = (-mrSges(7,1) * t156 - mrSges(7,2) * t159) * qJD(6) * pkin(5);
t116 = pkin(5) * t163 + t147;
t111 = Ifges(6,1) * t139 - Ifges(6,4) * t163;
t110 = Ifges(6,4) * t139 - Ifges(6,2) * t163;
t108 = mrSges(5,1) * t138 - mrSges(5,3) * t178;
t107 = -mrSges(5,2) * t138 - mrSges(5,3) * t179;
t103 = -Ifges(6,1) * t127 - Ifges(6,4) * t128;
t102 = -Ifges(6,4) * t127 - Ifges(6,2) * t128;
t101 = -t121 + t186;
t100 = mrSges(5,1) * t130 + mrSges(5,3) * t180;
t99 = -mrSges(5,2) * t130 + mrSges(5,3) * t181;
t82 = mrSges(6,1) * t138 + mrSges(6,3) * t96;
t81 = -mrSges(6,2) * t138 - mrSges(6,3) * t95;
t80 = t185 - (t154 * Ifges(5,1) - t189) * t129;
t79 = t184 - (-t183 + t188) * t129;
t73 = Ifges(7,1) * t106 + Ifges(7,4) * t105;
t72 = Ifges(7,4) * t106 + Ifges(7,2) * t105;
t71 = -mrSges(7,1) * t105 + mrSges(7,2) * t106;
t65 = pkin(5) * t95 + t91;
t53 = -Ifges(6,1) * t96 - Ifges(6,4) * t95 + Ifges(6,5) * t138;
t52 = -Ifges(6,4) * t96 - Ifges(6,2) * t95 + Ifges(6,6) * t138;
t51 = mrSges(7,1) * t138 - mrSges(7,3) * t64;
t50 = -mrSges(7,2) * t138 + mrSges(7,3) * t63;
t46 = -mrSges(6,2) * t130 + mrSges(6,3) * t62;
t45 = mrSges(6,1) * t130 - mrSges(6,3) * t61;
t36 = Ifges(7,1) * t69 + Ifges(7,4) * t70;
t35 = Ifges(7,4) * t69 + Ifges(7,2) * t70;
t33 = -t62 * pkin(5) + t76;
t32 = -mrSges(7,1) * t63 + mrSges(7,2) * t64;
t30 = Ifges(7,1) * t64 + Ifges(7,4) * t63 + Ifges(7,5) * t138;
t29 = Ifges(7,4) * t64 + Ifges(7,2) * t63 + Ifges(7,6) * t138;
t28 = Ifges(6,1) * t61 + Ifges(6,4) * t62 + t130 * Ifges(6,5);
t27 = Ifges(6,4) * t61 + Ifges(6,2) * t62 + t130 * Ifges(6,6);
t14 = -mrSges(7,2) * t130 + mrSges(7,3) * t19;
t13 = mrSges(7,1) * t130 - mrSges(7,3) * t18;
t7 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t130 * Ifges(7,5);
t6 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t130 * Ifges(7,6);
t1 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t153 ^ 2 + t155 ^ 2) * qJD(2) + (-0.2e1 * t89 * mrSges(4,3) - 0.2e1 * (-Ifges(4,4) + t165) * t129 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3) + Ifges(7,3)) * t130 + t170 + t201) * t138 + (-t152 * t79 + t154 * t80 + (-0.2e1 * Ifges(4,4) + t165) * t130 - (Ifges(5,1) * t150 + (2 * Ifges(4,1)) + (t183 - 0.2e1 * t188) * t152) * t129 + 0.2e1 * (mrSges(5,1) * t152 + mrSges(5,2) * t154 + mrSges(4,3)) * t90) * t140 + t130 * (-Ifges(6,5) * t96 - Ifges(6,6) * t95) + 0.2e1 * t76 * (mrSges(6,1) * t95 - mrSges(6,2) * t96) + 0.2e1 * t169 * t168 + (t47 * t74 + t48 * t75 - t187) * t200 + 0.2e1 * m(4) * (t115 * t89 - t187) + (t10 * t2 + t3 * t9 + t33 * t65) * t198 + (t11 * t26 + t12 * t25 + t76 * t91) * t199 + t130 * (Ifges(7,5) * t64 + Ifges(7,6) * t63) + 0.2e1 * (-t115 * t130 + t129 * t203) * mrSges(4,3) - 0.2e1 * t203 * t94 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + t19 * t29 + t18 * t30 + 0.2e1 * t33 * t32 + 0.2e1 * t25 * t45 + 0.2e1 * t26 * t46 + 0.2e1 * t2 * t50 + 0.2e1 * t3 * t51 + t61 * t53 + t62 * t52 + t63 * t6 + t64 * t7 + 0.2e1 * t65 * t8 + 0.2e1 * t11 * t81 + 0.2e1 * t12 * t82 + 0.2e1 * t91 * t31 - t95 * t27 - t96 * t28 + 0.2e1 * t75 * t99 + 0.2e1 * t74 * t100 + 0.2e1 * t48 * t107 + 0.2e1 * t47 * t108; t154 * t100 + t105 * t13 + t106 * t14 - t127 * t81 - t128 * t82 - t163 * t45 + t139 * t46 + t152 * t99 + t69 * t50 + t70 * t51 + m(7) * (t10 * t69 + t105 * t3 + t106 * t2 + t70 * t9) + m(6) * (t11 * t139 - t12 * t163 - t127 * t26 - t128 * t25) + m(5) * (t152 * t48 + t154 * t47) + t168; 0.2e1 * m(6) * (-t127 * t139 + t128 * t163) + 0.2e1 * m(7) * (t105 * t70 + t106 * t69); m(7) * (t10 * t23 + t116 * t33 + t65 * t193 + t2 * t56 + t24 * t9 + t3 * t55) + (t190 + t175) * t138 / 0.2e1 + (-t11 * t163 - t12 * t139 + t127 * t25 - t128 * t26) * mrSges(6,3) + (Ifges(6,5) * t139 + Ifges(7,5) * t106 - Ifges(6,6) * t163 + Ifges(7,6) * t105) * t130 / 0.2e1 + t76 * (mrSges(6,1) * t163 + mrSges(6,2) * t139) - t163 * t27 / 0.2e1 + (t10 * t70 + t2 * t105 - t3 * t106 - t9 * t69) * mrSges(7,3) + m(6) * (t11 * t114 + t112 * t12 + t147 * t76 + t25 * t88 + t26 * t87) - (t154 * (Ifges(5,1) * t152 + t188) / 0.2e1 - t152 * (Ifges(5,2) * t154 + t189) / 0.2e1 + Ifges(4,5)) * t129 - (-pkin(5) * t32 + t52 / 0.2e1) * t128 + m(5) * (-pkin(3) * t90 + (-t152 * t74 + t154 * t75) * qJD(4) + (-t47 * t152 + t48 * t154) * qJ(4)) + (qJ(4) * t99 + qJD(4) * t107 + t48 * mrSges(5,3) + t79 / 0.2e1 - t90 * mrSges(5,1) + t184 / 0.2e1) * t154 + (-qJ(4) * t100 - qJD(4) * t108 - t47 * mrSges(5,3) + t80 / 0.2e1 + t90 * mrSges(5,2) + t185 / 0.2e1) * t152 + t23 * t50 + t24 * t51 + t55 * t13 + t56 * t14 + t63 * t35 / 0.2e1 + t64 * t36 / 0.2e1 + t65 * t34 + t69 * t30 / 0.2e1 + t70 * t29 / 0.2e1 + t33 * t71 + t19 * t72 / 0.2e1 + t18 * t73 / 0.2e1 + t87 * t81 + t88 * t82 - t89 * mrSges(4,2) - t90 * mrSges(4,1) - pkin(3) * t94 + t91 * t101 - t95 * t102 / 0.2e1 - t96 * t103 / 0.2e1 + t105 * t6 / 0.2e1 + t106 * t7 / 0.2e1 + t62 * t110 / 0.2e1 + t61 * t111 / 0.2e1 + t112 * t45 + t114 * t46 + t116 * t8 - t127 * t53 / 0.2e1 - Ifges(4,6) * t130 + t139 * t28 / 0.2e1 + t147 * t31; m(6) * (-t112 * t128 - t114 * t127 + t139 * t87 - t163 * t88) + m(7) * (t105 * t24 + t106 * t23 + t55 * t70 + t56 * t69); 0.2e1 * t147 * t101 - t163 * t102 + t139 * t103 + t105 * t35 + t106 * t36 - t127 * t111 + 0.2e1 * t116 * t34 + t69 * t73 + t70 * t72 - (-0.2e1 * pkin(5) * t71 + t110) * t128 + (t116 * t193 + t23 * t56 + t24 * t55) * t198 + (t112 * t88 + t114 * t87) * t199 + 0.2e1 * (t105 * t23 - t106 * t24 - t55 * t69 + t56 * t70) * mrSges(7,3) + 0.2e1 * (t112 * t127 - t114 * t128 - t139 * t88 - t163 * t87) * mrSges(6,3) + (qJ(4) * t200 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t152 ^ 2 + t150); m(5) * t90 + m(6) * t76 + m(7) * t33 + t31 + t8 + t94; 0; -(-mrSges(6,1) - t196) * t128 + t166; 0; t12 * mrSges(6,1) - t11 * mrSges(6,2) + (m(7) * (t10 * t171 + t156 * t2 + t159 * t3 - t9 * t172) + t50 * t171 + t156 * t14 - t51 * t172 + t159 * t13) * pkin(5) + t162 + t201; -t186 + (t156 * t69 + t159 * t70 + (-t105 * t156 + t106 * t159) * qJD(6)) * t196 - t166; t88 * mrSges(6,1) - t87 * mrSges(6,2) + (m(7) * (t156 * t23 + t159 * t24 + (-t156 * t55 + t159 * t56) * qJD(6)) + (t156 * t70 - t159 * t69 + (t105 * t159 + t106 * t156) * qJD(6)) * mrSges(7,3)) * pkin(5) + t164 + t175; 0; 0.2e1 * t132; t162; -t34; t164; 0; t132; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
