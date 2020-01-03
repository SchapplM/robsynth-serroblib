% Calculate time derivative of joint inertia matrix for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:10
% EndTime: 2019-12-31 20:32:16
% DurationCPUTime: 2.11s
% Computational Cost: add. (3487->328), mult. (8351->512), div. (0->0), fcn. (7617->8), ass. (0->147)
t139 = sin(pkin(9));
t174 = pkin(7) + qJ(3);
t123 = t174 * t139;
t140 = cos(pkin(9));
t124 = t174 * t140;
t142 = sin(qJ(4));
t145 = cos(qJ(4));
t88 = -t142 * t123 + t145 * t124;
t148 = t139 * t142 - t140 * t145;
t105 = t148 * qJD(4);
t143 = sin(qJ(2));
t162 = qJD(2) * t143;
t118 = t139 * t145 + t140 * t142;
t106 = t118 * qJD(4);
t146 = cos(qJ(2));
t161 = qJD(2) * t146;
t67 = -t143 * t106 - t148 * t161;
t68 = t105 * t143 - t118 * t161;
t187 = -Ifges(5,5) * t67 - Ifges(5,6) * t68 - Ifges(5,3) * t162;
t186 = 2 * m(4);
t185 = 2 * m(5);
t184 = 2 * m(6);
t183 = -2 * pkin(1);
t182 = 2 * pkin(6);
t141 = sin(qJ(5));
t144 = cos(qJ(5));
t77 = -t118 * t141 - t144 * t148;
t181 = t77 / 0.2e1;
t78 = t118 * t144 - t141 * t148;
t180 = t78 / 0.2e1;
t179 = -t148 / 0.2e1;
t178 = t118 / 0.2e1;
t177 = t140 / 0.2e1;
t175 = pkin(4) * t106;
t39 = t77 * qJD(5) - t105 * t144 - t106 * t141;
t40 = -t78 * qJD(5) + t105 * t141 - t106 * t144;
t173 = Ifges(6,5) * t39 + Ifges(6,6) * t40;
t122 = -pkin(2) * t146 - t143 * qJ(3) - pkin(1);
t113 = t140 * t122;
t166 = t140 * t143;
t79 = -pkin(7) * t166 + t113 + (-pkin(6) * t139 - pkin(3)) * t146;
t168 = t139 * t143;
t165 = t140 * t146;
t94 = pkin(6) * t165 + t139 * t122;
t86 = -pkin(7) * t168 + t94;
t45 = t142 * t79 + t145 * t86;
t172 = mrSges(4,2) * t140;
t171 = Ifges(4,4) * t139;
t170 = Ifges(4,4) * t140;
t169 = t142 * t86;
t104 = -t143 * qJD(3) + (pkin(2) * t143 - qJ(3) * t146) * qJD(2);
t154 = pkin(6) * t162;
t84 = t140 * t104 + t139 * t154;
t167 = t139 * t146;
t164 = -Ifges(5,5) * t105 - Ifges(5,6) * t106;
t153 = t139 * t161;
t100 = mrSges(4,1) * t153 + t161 * t172;
t135 = pkin(6) * t161;
t109 = pkin(3) * t153 + t135;
t121 = pkin(3) * t168 + t143 * pkin(6);
t160 = qJD(3) * t139;
t159 = qJD(3) * t140;
t158 = qJD(4) * t145;
t157 = qJD(5) * t141;
t156 = qJD(5) * t144;
t98 = t118 * t143;
t99 = t148 * t143;
t59 = t141 * t99 - t144 * t98;
t24 = t59 * qJD(5) + t141 * t68 + t144 * t67;
t60 = -t141 * t98 - t144 * t99;
t25 = -t60 * qJD(5) - t141 * t67 + t144 * t68;
t155 = -Ifges(6,5) * t24 - Ifges(6,6) * t25 - Ifges(6,3) * t162;
t132 = -pkin(3) * t140 - pkin(2);
t35 = -t68 * mrSges(5,1) + t67 * mrSges(5,2);
t6 = -t25 * mrSges(6,1) + t24 * mrSges(6,2);
t16 = -mrSges(6,1) * t40 + t39 * mrSges(6,2);
t44 = t145 * t79 - t169;
t87 = -t145 * t123 - t124 * t142;
t65 = -pkin(8) * t118 + t87;
t66 = -pkin(8) * t148 + t88;
t32 = -t141 * t66 + t144 * t65;
t55 = -t123 * t158 + t145 * t159 + (-qJD(4) * t124 - t160) * t142;
t46 = -pkin(8) * t106 + t55;
t56 = -t118 * qJD(3) - qJD(4) * t88;
t47 = pkin(8) * t105 + t56;
t8 = t32 * qJD(5) + t141 * t47 + t144 * t46;
t33 = t141 * t65 + t144 * t66;
t9 = -t33 * qJD(5) - t141 * t46 + t144 * t47;
t152 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + t173;
t151 = t140 * Ifges(4,1) - t171;
t150 = -t139 * Ifges(4,2) + t170;
t149 = -Ifges(4,5) * t140 + Ifges(4,6) * t139;
t28 = -pkin(4) * t146 + t99 * pkin(8) + t44;
t34 = -pkin(8) * t98 + t45;
t12 = -t141 * t34 + t144 * t28;
t13 = t141 * t28 + t144 * t34;
t61 = (pkin(3) * t143 - pkin(7) * t165) * qJD(2) + t84;
t96 = t139 * t104;
t70 = t96 + (-pkin(6) * t166 - pkin(7) * t167) * qJD(2);
t15 = -t45 * qJD(4) - t142 * t70 + t145 * t61;
t10 = pkin(4) * t162 - pkin(8) * t67 + t15;
t14 = -qJD(4) * t169 + t142 * t61 + t145 * t70 + t79 * t158;
t11 = pkin(8) * t68 + t14;
t2 = t12 * qJD(5) + t10 * t141 + t11 * t144;
t3 = -t13 * qJD(5) + t10 * t144 - t11 * t141;
t147 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t155;
t120 = -mrSges(4,1) * t146 - mrSges(4,3) * t166;
t119 = mrSges(4,2) * t146 - mrSges(4,3) * t168;
t111 = (-mrSges(6,1) * t141 - mrSges(6,2) * t144) * qJD(5) * pkin(4);
t108 = (mrSges(4,1) * t143 - mrSges(4,3) * t165) * qJD(2);
t107 = (-mrSges(4,2) * t143 - mrSges(4,3) * t167) * qJD(2);
t101 = t105 * mrSges(5,2);
t95 = pkin(4) * t148 + t132;
t93 = -pkin(6) * t167 + t113;
t92 = (t143 * Ifges(4,5) + t151 * t146) * qJD(2);
t91 = (t143 * Ifges(4,6) + t150 * t146) * qJD(2);
t90 = -mrSges(5,1) * t146 + t99 * mrSges(5,3);
t89 = mrSges(5,2) * t146 - t98 * mrSges(5,3);
t85 = -t140 * t154 + t96;
t83 = Ifges(5,1) * t118 - Ifges(5,4) * t148;
t82 = Ifges(5,4) * t118 - Ifges(5,2) * t148;
t81 = pkin(4) * t98 + t121;
t76 = -Ifges(5,1) * t105 - Ifges(5,4) * t106;
t75 = -Ifges(5,4) * t105 - Ifges(5,2) * t106;
t74 = mrSges(5,1) * t106 - t101;
t58 = -Ifges(5,1) * t99 - Ifges(5,4) * t98 - Ifges(5,5) * t146;
t57 = -Ifges(5,4) * t99 - Ifges(5,2) * t98 - Ifges(5,6) * t146;
t52 = -mrSges(5,2) * t162 + mrSges(5,3) * t68;
t51 = mrSges(5,1) * t162 - mrSges(5,3) * t67;
t50 = -mrSges(6,1) * t146 - t60 * mrSges(6,3);
t49 = mrSges(6,2) * t146 + t59 * mrSges(6,3);
t48 = -pkin(4) * t68 + t109;
t43 = Ifges(6,1) * t78 + Ifges(6,4) * t77;
t42 = Ifges(6,4) * t78 + Ifges(6,2) * t77;
t41 = -mrSges(6,1) * t77 + mrSges(6,2) * t78;
t31 = -mrSges(6,1) * t59 + mrSges(6,2) * t60;
t30 = Ifges(5,1) * t67 + Ifges(5,4) * t68 + Ifges(5,5) * t162;
t29 = Ifges(5,4) * t67 + Ifges(5,2) * t68 + Ifges(5,6) * t162;
t27 = Ifges(6,1) * t60 + Ifges(6,4) * t59 - Ifges(6,5) * t146;
t26 = Ifges(6,4) * t60 + Ifges(6,2) * t59 - Ifges(6,6) * t146;
t20 = -mrSges(6,2) * t162 + mrSges(6,3) * t25;
t19 = mrSges(6,1) * t162 - mrSges(6,3) * t24;
t18 = Ifges(6,1) * t39 + Ifges(6,4) * t40;
t17 = Ifges(6,4) * t39 + Ifges(6,2) * t40;
t5 = Ifges(6,1) * t24 + Ifges(6,4) * t25 + Ifges(6,5) * t162;
t4 = Ifges(6,4) * t24 + Ifges(6,2) * t25 + Ifges(6,6) * t162;
t1 = [(((mrSges(3,2) * t183) + 0.2e1 * (Ifges(3,4) + t149) * t146) * t146 + (-Ifges(5,5) * t99 - Ifges(5,6) * t98 + Ifges(6,5) * t60 + Ifges(6,6) * t59 + (mrSges(3,1) * t183) + (-0.2e1 * Ifges(3,4) - t149) * t143 + ((mrSges(4,1) * t139 + t172) * t182 - (2 * Ifges(4,3)) + (pkin(6) ^ 2 * t186) - Ifges(5,3) - Ifges(6,3) - (2 * Ifges(3,2)) + (2 * Ifges(3,1)) + t140 * t151 - t139 * t150) * t146) * t143) * qJD(2) + 0.2e1 * t109 * (t98 * mrSges(5,1) - t99 * mrSges(5,2)) + (t155 + t187) * t146 + (t12 * t3 + t13 * t2 + t48 * t81) * t184 + (t109 * t121 + t14 * t45 + t15 * t44) * t185 + (t93 * t84 + t94 * t85) * t186 + (t100 * t182 - t139 * t91 + t140 * t92) * t143 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t20 + t25 * t26 + t24 * t27 + 0.2e1 * t48 * t31 + 0.2e1 * t2 * t49 + 0.2e1 * t3 * t50 + 0.2e1 * t44 * t51 + 0.2e1 * t45 * t52 + t59 * t4 + t60 * t5 + t67 * t58 + t68 * t57 + 0.2e1 * t81 * t6 + 0.2e1 * t14 * t89 + 0.2e1 * t15 * t90 - t98 * t29 - t99 * t30 + 0.2e1 * t94 * t107 + 0.2e1 * t93 * t108 + 0.2e1 * t85 * t119 + 0.2e1 * t84 * t120 + 0.2e1 * t121 * t35; (qJD(3) * t119 + qJ(3) * t107 + t85 * mrSges(4,3) + t91 / 0.2e1) * t140 + (-qJD(3) * t120 - qJ(3) * t108 - t84 * mrSges(4,3) + t92 / 0.2e1) * t139 - (-pkin(4) * t31 + t57 / 0.2e1) * t106 + (-t12 * t39 + t13 * t40 + t2 * t77 - t3 * t78) * mrSges(6,3) - (t173 + t164) * t146 / 0.2e1 + m(6) * (t12 * t9 + t13 * t8 + t81 * t175 + t2 * t33 + t3 * t32 + t48 * t95) + m(4) * (t94 * t159 - t93 * t160 + (-t84 * t139 + t85 * t140) * qJ(3)) + (t44 * t105 - t45 * t106 - t15 * t118 - t14 * t148) * mrSges(5,3) + t109 * (mrSges(5,1) * t148 + t118 * mrSges(5,2)) + t30 * t178 + t29 * t179 + t5 * t180 + t4 * t181 + m(5) * (t109 * t132 + t14 * t88 + t15 * t87 + t44 * t56 + t45 * t55) + (((pkin(6) * mrSges(3,2)) + Ifges(4,5) * t139 / 0.2e1 + Ifges(4,6) * t177 + Ifges(5,5) * t178 + Ifges(5,6) * t179 + Ifges(6,5) * t180 + Ifges(6,6) * t181 - Ifges(3,6)) * t143 + ((Ifges(4,1) * t139 + t170) * t177 - t139 * (Ifges(4,2) * t140 + t171) / 0.2e1 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(4,1) * t140 + mrSges(4,2) * t139 - mrSges(3,1)) * pkin(6)) * t146) * qJD(2) + t32 * t19 + t33 * t20 + t39 * t27 / 0.2e1 + t40 * t26 / 0.2e1 + t25 * t42 / 0.2e1 + t24 * t43 / 0.2e1 + t48 * t41 + t8 * t49 + t9 * t50 + t59 * t17 / 0.2e1 + t60 * t18 / 0.2e1 + t81 * t16 + t68 * t82 / 0.2e1 + t67 * t83 / 0.2e1 + t87 * t51 + t88 * t52 + t55 * t89 + t56 * t90 + t95 * t6 - t98 * t75 / 0.2e1 - t99 * t76 / 0.2e1 - pkin(2) * t100 - t105 * t58 / 0.2e1 + t121 * t74 + t132 * t35; -t105 * t83 - t148 * t75 + t118 * t76 + 0.2e1 * t132 * t74 + 0.2e1 * t95 * t16 + t77 * t17 + t78 * t18 + t39 * t43 + t40 * t42 - (-0.2e1 * pkin(4) * t41 + t82) * t106 + (t95 * t175 + t32 * t9 + t33 * t8) * t184 + (t55 * t88 + t56 * t87) * t185 + 0.2e1 * (-t32 * t39 + t33 * t40 + t77 * t8 - t78 * t9) * mrSges(6,3) + 0.2e1 * (t105 * t87 - t106 * t88 - t118 * t56 - t148 * t55) * mrSges(5,3) + (qJ(3) * t186 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t139 ^ 2 + t140 ^ 2); m(4) * t135 + m(5) * t109 + m(6) * t48 + t100 + t35 + t6; -t101 - (-m(6) * pkin(4) - mrSges(5,1)) * t106 + t16; 0; t15 * mrSges(5,1) - t14 * mrSges(5,2) + (m(6) * (-t12 * t157 + t13 * t156 + t141 * t2 + t144 * t3) + t49 * t156 + t141 * t20 - t50 * t157 + t144 * t19) * pkin(4) + t147 - t187; t56 * mrSges(5,1) - t55 * mrSges(5,2) + (m(6) * (t141 * t8 + t144 * t9 + (-t141 * t32 + t144 * t33) * qJD(5)) + (t141 * t40 - t144 * t39 + (t141 * t78 + t144 * t77) * qJD(5)) * mrSges(6,3)) * pkin(4) + t152 + t164; 0; 0.2e1 * t111; t147; t152; 0; t111; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
