% Calculate time derivative of joint inertia matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:00
% EndTime: 2019-12-31 21:00:08
% DurationCPUTime: 2.94s
% Computational Cost: add. (1832->336), mult. (4721->488), div. (0->0), fcn. (3695->6), ass. (0->148)
t193 = Ifges(6,4) + Ifges(5,5);
t192 = Ifges(6,6) - Ifges(5,6);
t130 = sin(qJ(3));
t131 = sin(qJ(2));
t132 = cos(qJ(3));
t157 = qJD(3) * t132;
t133 = cos(qJ(2));
t160 = qJD(2) * t133;
t191 = -t130 * t160 - t131 * t157;
t195 = -Ifges(6,2) - Ifges(4,3) - Ifges(5,3);
t129 = sin(pkin(8));
t167 = cos(pkin(8));
t146 = t167 * t130;
t95 = t129 * t132 + t146;
t89 = t95 * qJD(3);
t145 = t167 * t132;
t166 = t129 * t130;
t139 = t145 - t166;
t90 = t139 * qJD(3);
t194 = Ifges(4,5) * t157 + t192 * t89 + t193 * t90;
t150 = t132 * t160;
t158 = qJD(3) * t131;
t152 = t130 * t158;
t137 = t150 - t152;
t180 = pkin(3) * t129;
t119 = qJ(5) + t180;
t190 = m(6) * t119 + mrSges(6,3);
t189 = 2 * m(4);
t188 = 2 * m(5);
t187 = 0.2e1 * m(6);
t186 = -2 * pkin(1);
t185 = 2 * pkin(6);
t184 = m(5) * pkin(3);
t183 = -t130 / 0.2e1;
t182 = t132 / 0.2e1;
t179 = pkin(6) * t130;
t178 = t90 * mrSges(6,2);
t176 = -qJ(4) - pkin(7);
t109 = -pkin(2) * t133 - t131 * pkin(7) - pkin(1);
t162 = t132 * t133;
t118 = pkin(6) * t162;
t156 = qJD(4) * t132;
t104 = (pkin(2) * t131 - pkin(7) * t133) * qJD(2);
t161 = qJD(2) * t131;
t168 = t132 * t104 + t161 * t179;
t13 = -t131 * t156 + (pkin(3) * t131 - qJ(4) * t162) * qJD(2) + (-t118 + (qJ(4) * t131 - t109) * t130) * qJD(3) + t168;
t164 = t131 * t132;
t173 = t130 * t104 + t109 * t157;
t17 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t164 + (-qJD(4) * t131 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t133) * t130 + t173;
t4 = t129 * t13 + t167 * t17;
t97 = t132 * t109;
t57 = -qJ(4) * t164 + t97 + (-pkin(3) - t179) * t133;
t165 = t130 * t131;
t76 = t130 * t109 + t118;
t65 = -qJ(4) * t165 + t76;
t20 = t129 * t57 + t167 * t65;
t172 = Ifges(4,4) * t130;
t171 = Ifges(4,4) * t132;
t170 = Ifges(4,6) * t133;
t143 = Ifges(4,1) * t132 - t172;
t79 = -Ifges(4,5) * t133 + t131 * t143;
t169 = t132 * t79;
t112 = Ifges(4,1) * t130 + t171;
t163 = t132 * t112;
t105 = pkin(3) * t165 + t131 * pkin(6);
t159 = qJD(3) * t130;
t155 = pkin(3) * t159;
t74 = -t191 * pkin(3) + pkin(6) * t160;
t122 = -pkin(3) * t132 - pkin(2);
t154 = t167 * pkin(3);
t40 = -t129 * t137 - t145 * t158 - t146 * t160;
t41 = t139 * t160 - t158 * t95;
t11 = -t40 * mrSges(5,1) + t41 * mrSges(5,2);
t49 = t89 * mrSges(5,1) + t90 * mrSges(5,2);
t10 = -t40 * mrSges(6,1) - t41 * mrSges(6,3);
t48 = t89 * mrSges(6,1) - t90 * mrSges(6,3);
t147 = qJD(3) * t176;
t135 = -qJD(4) * t130 + t132 * t147;
t88 = t130 * t147 + t156;
t42 = t129 * t88 - t135 * t167;
t43 = t129 * t135 + t167 * t88;
t110 = t176 * t132;
t66 = -t110 * t129 - t146 * t176;
t67 = -t110 * t167 + t166 * t176;
t149 = t42 * t66 + t67 * t43;
t148 = Ifges(4,6) * t130 + (2 * Ifges(3,4));
t26 = -mrSges(6,1) * t161 + t41 * mrSges(6,2);
t144 = mrSges(4,1) * t130 + mrSges(4,2) * t132;
t142 = -Ifges(4,2) * t130 + t171;
t111 = Ifges(4,2) * t132 + t172;
t141 = Ifges(4,5) * t130 + Ifges(4,6) * t132;
t27 = (-t132 * t161 - t133 * t159) * pkin(6) + t173;
t28 = -t76 * qJD(3) + t168;
t140 = -t130 * t28 + t132 * t27;
t3 = -t129 * t17 + t13 * t167;
t19 = -t129 * t65 + t167 * t57;
t138 = -Ifges(4,5) * t150 + t195 * t161 + t192 * t40 - t193 * t41;
t121 = -t154 - pkin(4);
t103 = -mrSges(4,1) * t133 - mrSges(4,3) * t164;
t102 = mrSges(4,2) * t133 - mrSges(4,3) * t165;
t101 = t143 * qJD(3);
t100 = t142 * qJD(3);
t99 = t144 * qJD(3);
t81 = t139 * t131;
t80 = t95 * t131;
t78 = t131 * t142 - t170;
t75 = -t133 * t179 + t97;
t73 = -mrSges(4,2) * t161 + mrSges(4,3) * t191;
t72 = mrSges(4,1) * t161 - mrSges(4,3) * t137;
t71 = mrSges(6,1) * t133 + t81 * mrSges(6,2);
t70 = -mrSges(5,1) * t133 - t81 * mrSges(5,3);
t69 = mrSges(5,2) * t133 - t80 * mrSges(5,3);
t68 = -t80 * mrSges(6,2) - mrSges(6,3) * t133;
t64 = Ifges(5,1) * t95 + Ifges(5,4) * t139;
t63 = Ifges(6,1) * t95 - Ifges(6,5) * t139;
t62 = Ifges(5,4) * t95 + Ifges(5,2) * t139;
t61 = Ifges(6,5) * t95 - Ifges(6,3) * t139;
t60 = -mrSges(5,1) * t139 + mrSges(5,2) * t95;
t59 = -mrSges(6,1) * t139 - mrSges(6,3) * t95;
t56 = -mrSges(4,1) * t191 + mrSges(4,2) * t137;
t54 = -pkin(4) * t139 - qJ(5) * t95 + t122;
t53 = Ifges(5,1) * t90 - Ifges(5,4) * t89;
t52 = Ifges(6,1) * t90 + Ifges(6,5) * t89;
t51 = Ifges(5,4) * t90 - Ifges(5,2) * t89;
t50 = Ifges(6,5) * t90 + Ifges(6,3) * t89;
t47 = -t112 * t158 + (Ifges(4,5) * t131 + t133 * t143) * qJD(2);
t46 = -t111 * t158 + (Ifges(4,6) * t131 + t133 * t142) * qJD(2);
t45 = mrSges(5,1) * t80 + mrSges(5,2) * t81;
t44 = mrSges(6,1) * t80 - mrSges(6,3) * t81;
t32 = Ifges(5,1) * t81 - Ifges(5,4) * t80 - Ifges(5,5) * t133;
t31 = Ifges(6,1) * t81 - Ifges(6,4) * t133 + Ifges(6,5) * t80;
t30 = Ifges(5,4) * t81 - Ifges(5,2) * t80 - Ifges(5,6) * t133;
t29 = Ifges(6,5) * t81 - Ifges(6,6) * t133 + Ifges(6,3) * t80;
t25 = mrSges(5,1) * t161 - mrSges(5,3) * t41;
t24 = -mrSges(5,2) * t161 + mrSges(5,3) * t40;
t23 = mrSges(6,2) * t40 + mrSges(6,3) * t161;
t22 = pkin(4) * t80 - qJ(5) * t81 + t105;
t21 = pkin(4) * t89 - qJ(5) * t90 - qJD(5) * t95 + t155;
t16 = t133 * pkin(4) - t19;
t15 = -qJ(5) * t133 + t20;
t9 = Ifges(5,1) * t41 + Ifges(5,4) * t40 + Ifges(5,5) * t161;
t8 = Ifges(6,1) * t41 + Ifges(6,4) * t161 - Ifges(6,5) * t40;
t7 = Ifges(5,4) * t41 + Ifges(5,2) * t40 + Ifges(5,6) * t161;
t6 = Ifges(6,5) * t41 + Ifges(6,6) * t161 - Ifges(6,3) * t40;
t5 = -pkin(4) * t40 - qJ(5) * t41 - qJD(5) * t81 + t74;
t2 = -pkin(4) * t161 - t3;
t1 = qJ(5) * t161 - qJD(5) * t133 + t4;
t12 = [(t56 * t185 - t130 * t46 + t132 * t47 + (-t130 * t79 - t132 * t78 + t133 * t141) * qJD(3) + ((mrSges(3,1) * t186) + (Ifges(4,5) * t132 - t148) * t131 + t193 * t81 + t192 * t80 + ((pkin(6) ^ 2 * t189) + t144 * t185 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + t195) * t133) * qJD(2)) * t131 + (t8 + t9) * t81 + (t6 - t7) * t80 + (((mrSges(3,2) * t186) - t130 * t78 + t133 * t148 + t169) * qJD(2) + t138) * t133 + 0.2e1 * t27 * t102 + 0.2e1 * t28 * t103 + 0.2e1 * t105 * t11 + 0.2e1 * t75 * t72 + 0.2e1 * t76 * t73 + 0.2e1 * t1 * t68 + 0.2e1 * t4 * t69 + 0.2e1 * t3 * t70 + 0.2e1 * t2 * t71 + 0.2e1 * t74 * t45 + 0.2e1 * t5 * t44 + 0.2e1 * t19 * t25 + 0.2e1 * t16 * t26 + 0.2e1 * t22 * t10 + 0.2e1 * t15 * t23 + 0.2e1 * t20 * t24 + (t31 + t32) * t41 + (t30 - t29) * t40 + (t1 * t15 + t16 * t2 + t22 * t5) * t187 + (t105 * t74 + t19 * t3 + t20 * t4) * t188 + (t76 * t27 + t75 * t28) * t189; -t194 * t133 / 0.2e1 + (t23 + t24) * t67 + (t26 - t25) * t66 + (t68 + t69) * t43 + (-t70 + t71) * t42 + (Ifges(3,5) + t163 / 0.2e1 + t111 * t183 + (-m(4) * pkin(2) - mrSges(4,1) * t132 + mrSges(4,2) * t130 - mrSges(3,1)) * pkin(6)) * t160 + t130 * t47 / 0.2e1 + t122 * t11 + m(6) * (t1 * t67 + t15 * t43 + t16 * t42 + t2 * t66 + t21 * t22 + t5 * t54) + (t31 / 0.2e1 + t32 / 0.2e1) * t90 + (t29 / 0.2e1 - t30 / 0.2e1) * t89 + (t52 / 0.2e1 + t53 / 0.2e1) * t81 + (t50 / 0.2e1 - t51 / 0.2e1) * t80 + t105 * t49 + (t169 / 0.2e1 + (t170 / 0.2e1 - t78 / 0.2e1 + pkin(3) * t45) * t130) * qJD(3) + t5 * t59 + (t132 * t73 + m(4) * (-t157 * t75 - t159 * t76 + t140) - t130 * t72 - t103 * t157 - t102 * t159) * pkin(7) + t74 * t60 + t21 * t44 + t22 * t48 + t54 * t10 - pkin(2) * t56 + m(5) * (t105 * t155 + t122 * t74 - t19 * t42 + t20 * t43 - t3 * t66 + t4 * t67) + (t63 / 0.2e1 + t64 / 0.2e1) * t41 + (t8 / 0.2e1 + t9 / 0.2e1) * t95 + (-t61 / 0.2e1 + t62 / 0.2e1) * t40 + ((-t130 * t76 - t132 * t75) * qJD(3) + t140) * mrSges(4,3) - (t6 / 0.2e1 - t7 / 0.2e1) * t139 + (t139 * t4 - t19 * t90 - t20 * t89 - t3 * t95) * mrSges(5,3) + (t1 * t139 - t15 * t89 + t16 * t90 + t2 * t95) * mrSges(6,2) + (t101 * t182 + t100 * t183 - Ifges(3,6) * qJD(2) + (-t132 * t111 / 0.2e1 + t112 * t183) * qJD(3) + (qJD(2) * mrSges(3,2) + t99) * pkin(6) + (-t139 * t192 + t193 * t95 + t141) * qJD(2) / 0.2e1) * t131 + t46 * t182; -0.2e1 * pkin(2) * t99 + t132 * t100 + t130 * t101 + 0.2e1 * t122 * t49 + 0.2e1 * t21 * t59 + 0.2e1 * t54 * t48 + (t52 + t53) * t95 - (t50 - t51) * t139 + (t63 + t64) * t90 + (t61 - t62) * t89 + (t163 + (0.2e1 * pkin(3) * t60 - t111) * t130) * qJD(3) + (t122 * t155 + t149) * t188 + (t21 * t54 + t149) * t187 + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * (t139 * t43 + t42 * t95 + t66 * t90 - t67 * t89); t119 * t23 + t121 * t26 + t24 * t180 + qJD(5) * t68 - t27 * mrSges(4,2) + t28 * mrSges(4,1) + m(6) * (qJD(5) * t15 + t1 * t119 + t121 * t2) + t3 * mrSges(5,1) - t4 * mrSges(5,2) + t1 * mrSges(6,3) - t2 * mrSges(6,1) + t25 * t154 - Ifges(4,5) * t152 - t138 + (t129 * t4 + t167 * t3) * t184 + t191 * Ifges(4,6); m(6) * qJD(5) * t67 - Ifges(4,6) * t159 + t121 * t178 + (t129 * t184 - mrSges(5,2) + t190) * t43 + (m(6) * t121 - t167 * t184 - mrSges(5,1) - mrSges(6,1)) * t42 + (-mrSges(4,1) * t157 + mrSges(4,2) * t159) * pkin(7) + (-t154 * t90 - t180 * t89) * mrSges(5,3) + (qJD(5) * t139 - t119 * t89) * mrSges(6,2) + t194; 0.2e1 * t190 * qJD(5); m(5) * t74 + m(6) * t5 + t10 + t11; m(5) * t155 + m(6) * t21 + t48 + t49; 0; 0; m(6) * t2 + t26; m(6) * t42 + t178; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;
