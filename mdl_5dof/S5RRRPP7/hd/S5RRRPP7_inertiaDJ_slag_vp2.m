% Calculate time derivative of joint inertia matrix for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:58
% EndTime: 2019-12-31 21:04:03
% DurationCPUTime: 2.04s
% Computational Cost: add. (1033->335), mult. (2649->459), div. (0->0), fcn. (1687->4), ass. (0->149)
t173 = Ifges(5,4) + Ifges(4,5);
t174 = -Ifges(5,2) - Ifges(4,3);
t98 = cos(qJ(2));
t139 = qJD(2) * t98;
t97 = cos(qJ(3));
t123 = t97 * t139;
t96 = sin(qJ(2));
t137 = qJD(3) * t96;
t95 = sin(qJ(3));
t126 = t95 * t137;
t101 = t123 - t126;
t136 = qJD(3) * t97;
t100 = t96 * t136 + t95 * t139;
t154 = Ifges(5,5) * t97;
t110 = Ifges(5,3) * t95 + t154;
t156 = Ifges(6,4) * t97;
t111 = Ifges(6,2) * t95 + t156;
t172 = (t110 + t111) * qJD(3);
t140 = qJD(2) * t96;
t171 = qJ(4) * t140 - qJD(4) * t98;
t157 = Ifges(6,4) * t95;
t113 = Ifges(6,1) * t97 + t157;
t155 = Ifges(5,5) * t95;
t114 = Ifges(5,1) * t97 + t155;
t159 = Ifges(4,4) * t95;
t115 = Ifges(4,1) * t97 - t159;
t170 = (t113 + t114 + t115) * qJD(3);
t143 = qJ(4) * t95;
t161 = pkin(3) + pkin(4);
t169 = -t161 * t97 - t143;
t168 = 2 * m(4);
t167 = 2 * m(6);
t166 = -2 * pkin(1);
t165 = 2 * pkin(6);
t164 = -2 * mrSges(6,3);
t160 = pkin(6) * t95;
t158 = Ifges(4,4) * t97;
t153 = t95 * t96;
t152 = t96 * t97;
t151 = t97 * t98;
t150 = t98 * Ifges(6,5);
t148 = -mrSges(5,2) + mrSges(6,3);
t147 = Ifges(4,6) + Ifges(6,6);
t146 = pkin(7) - qJ(5);
t64 = (pkin(2) * t96 - pkin(7) * t98) * qJD(2);
t66 = -pkin(2) * t98 - pkin(7) * t96 - pkin(1);
t145 = t66 * t136 + t95 * t64;
t138 = qJD(3) * t95;
t86 = pkin(6) * t151;
t144 = qJD(3) * t86 + t66 * t138;
t31 = t95 * t66 + t86;
t142 = qJ(4) * t97;
t141 = qJ(5) * t96;
t135 = qJD(4) * t95;
t134 = qJD(4) * t97;
t132 = qJD(5) * t97;
t131 = qJ(5) * qJD(2);
t130 = qJ(5) * qJD(3);
t129 = mrSges(6,3) * t151;
t128 = -Ifges(6,5) + t173;
t35 = t113 * t96 + t150;
t36 = -Ifges(5,4) * t98 + t114 * t96;
t37 = -Ifges(4,5) * t98 + t115 * t96;
t127 = -t35 - t36 - t37;
t122 = -pkin(3) - t160;
t70 = t146 * t97;
t85 = t98 * t160;
t30 = t66 * t97 - t85;
t121 = -t64 * t97 + t144;
t21 = -qJ(4) * t98 + t31;
t71 = -Ifges(5,3) * t97 + t155;
t72 = -Ifges(6,2) * t97 + t157;
t73 = Ifges(4,2) * t97 + t159;
t120 = t71 / 0.2e1 + t72 / 0.2e1 - t73 / 0.2e1;
t74 = Ifges(6,1) * t95 - t156;
t75 = Ifges(5,1) * t95 - t154;
t76 = Ifges(4,1) * t95 + t158;
t119 = t74 / 0.2e1 + t75 / 0.2e1 + t76 / 0.2e1;
t50 = -mrSges(6,1) * t138 + mrSges(6,2) * t136;
t118 = -t97 * mrSges(4,1) + t95 * mrSges(4,2);
t117 = mrSges(4,1) * t95 + mrSges(4,2) * t97;
t68 = -t97 * mrSges(5,1) - t95 * mrSges(5,3);
t116 = mrSges(5,1) * t95 - mrSges(5,3) * t97;
t112 = -Ifges(4,2) * t95 + t158;
t109 = pkin(3) * t97 + t143;
t108 = pkin(3) * t95 - t142;
t107 = (m(5) * pkin(7) - t148) * t97;
t106 = -t100 * Ifges(5,6) - t173 * t123 + t174 * t140;
t105 = pkin(6) + t108;
t104 = -t161 * t95 + t142;
t32 = -Ifges(5,6) * t98 + t110 * t96;
t33 = t98 * Ifges(6,6) + t111 * t96;
t34 = -Ifges(4,6) * t98 + t112 * t96;
t103 = t147 * t98 + t32 + t33 - t34;
t102 = -pkin(6) + t104;
t7 = (-t138 * t98 - t140 * t97) * pkin(6) + t145;
t17 = -mrSges(6,1) * t100 + t101 * mrSges(6,2);
t94 = t98 * pkin(3);
t93 = Ifges(5,4) * t136;
t92 = Ifges(4,5) * t136;
t90 = Ifges(5,6) * t138;
t79 = mrSges(6,3) * t126;
t78 = mrSges(5,2) * t123;
t69 = mrSges(6,1) * t97 + mrSges(6,2) * t95;
t67 = t146 * t95;
t65 = -pkin(2) - t109;
t63 = -mrSges(5,2) * t153 - mrSges(5,3) * t98;
t62 = mrSges(5,1) * t98 + mrSges(5,2) * t152;
t61 = -mrSges(4,1) * t98 - mrSges(4,3) * t152;
t60 = mrSges(6,1) * t98 - mrSges(6,3) * t152;
t59 = mrSges(4,2) * t98 - mrSges(4,3) * t153;
t58 = -mrSges(6,2) * t98 + mrSges(6,3) * t153;
t54 = t112 * qJD(3);
t51 = t117 * qJD(3);
t49 = t116 * qJD(3);
t47 = pkin(2) - t169;
t44 = (-mrSges(6,1) * t95 + mrSges(6,2) * t97) * t96;
t43 = t116 * t96;
t41 = qJD(3) * t70 - qJD(5) * t95;
t40 = qJD(3) * t108 - t135;
t39 = -t138 * t146 - t132;
t38 = t105 * t96;
t29 = qJD(3) * t104 + t135;
t28 = -mrSges(5,2) * t100 + mrSges(5,3) * t140;
t27 = -mrSges(4,2) * t140 - mrSges(4,3) * t100;
t26 = mrSges(6,2) * t140 + mrSges(6,3) * t100;
t25 = t78 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t138) * t96;
t24 = mrSges(4,1) * t140 - mrSges(4,3) * t101;
t23 = t79 + (-mrSges(6,1) * t96 - t129) * qJD(2);
t22 = -t30 + t94;
t20 = t102 * t96;
t19 = t141 * t95 + t21;
t18 = mrSges(4,1) * t100 + mrSges(4,2) * t101;
t16 = mrSges(5,1) * t100 - mrSges(5,3) * t101;
t15 = pkin(4) * t98 + t85 + t94 + (-t66 - t141) * t97;
t14 = -t76 * t137 + (Ifges(4,5) * t96 + t115 * t98) * qJD(2);
t13 = -t75 * t137 + (Ifges(5,4) * t96 + t114 * t98) * qJD(2);
t12 = -t74 * t137 + (-Ifges(6,5) * t96 + t113 * t98) * qJD(2);
t11 = -t73 * t137 + (Ifges(4,6) * t96 + t112 * t98) * qJD(2);
t10 = -t72 * t137 + (-Ifges(6,6) * t96 + t111 * t98) * qJD(2);
t9 = -t71 * t137 + (Ifges(5,6) * t96 + t110 * t98) * qJD(2);
t8 = t140 * t160 - t121;
t6 = (qJD(3) * t109 - t134) * t96 + t105 * t139;
t5 = t122 * t140 + t121;
t4 = t7 + t171;
t3 = (t169 * qJD(3) + t134) * t96 + t102 * t139;
t2 = (-qJD(2) * pkin(6) + t130) * t152 + (qJD(5) * t96 + (-pkin(6) * qJD(3) + t131) * t98) * t95 + t145 + t171;
t1 = (-t131 * t98 - t64) * t97 + (t95 * t130 - t132 + (-pkin(4) + t122) * qJD(2)) * t96 + t144;
t42 = [0.2e1 * t2 * t58 + 0.2e1 * t7 * t59 + 0.2e1 * t1 * t60 + 0.2e1 * t8 * t61 + 0.2e1 * t5 * t62 + 0.2e1 * t4 * t63 + 0.2e1 * t31 * t27 + 0.2e1 * t38 * t16 + 0.2e1 * t6 * t43 + 0.2e1 * t3 * t44 + 0.2e1 * t20 * t17 + 0.2e1 * t15 * t23 + 0.2e1 * t22 * t25 + 0.2e1 * t19 * t26 + 0.2e1 * t21 * t28 + 0.2e1 * t30 * t24 + (t30 * t8 + t31 * t7) * t168 + 0.2e1 * m(5) * (t21 * t4 + t22 * t5 + t38 * t6) + (t1 * t15 + t19 * t2 + t20 * t3) * t167 + (((mrSges(3,2) * t166) + 0.2e1 * Ifges(3,4) * t98 + (-t127 + t150) * t97 + t103 * t95) * qJD(2) + t106) * t98 + (t18 * t165 + (t12 + t13 + t14) * t97 + (t10 - t11 + t9) * t95 + (t103 * t97 + (t128 * t98 + t127) * t95) * qJD(3) + ((mrSges(3,1) * t166) - 0.2e1 * Ifges(3,4) * t96 + t128 * t152 + (Ifges(5,6) - t147) * t153 + ((pkin(6) ^ 2 * t168) + t117 * t165 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(6,3)) + t174) * t98) * qJD(2)) * t96; -pkin(2) * t18 + t65 * t16 + t47 * t17 + t20 * t50 + t67 * t23 + t70 * t26 + t29 * t44 + t3 * t69 + t38 * t49 + t39 * t58 + t40 * t43 + t41 * t60 + t6 * t68 + m(6) * (t1 * t67 + t15 * t41 + t19 * t39 + t2 * t70 + t20 * t29 + t3 * t47) + m(5) * (t38 * t40 + t6 * t65) + (t7 * mrSges(4,3) - t2 * mrSges(6,3) + t4 * mrSges(5,2) - t9 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1) * t97 + (-t8 * mrSges(4,3) - t1 * mrSges(6,3) + t5 * mrSges(5,2) + t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1) * t95 + ((t27 + t28) * t97 + (-t24 + t25) * t95 + m(5) * (t4 * t97 + t5 * t95) + m(4) * (t7 * t97 - t8 * t95)) * pkin(7) + (-t92 / 0.2e1 - t93 / 0.2e1 - t90 / 0.2e1 + (Ifges(3,5) + t119 * t97 + t120 * t95 + (-m(4) * pkin(2) - mrSges(3,1) + t118) * pkin(6)) * qJD(2)) * t98 + (-qJD(2) * (Ifges(6,5) * t95 - Ifges(6,6) * t97) / 0.2e1 - Ifges(3,6) * qJD(2) - t95 * t54 / 0.2e1 + (qJD(2) * mrSges(3,2) + t51) * pkin(6) + t172 * t95 / 0.2e1 + ((Ifges(4,6) - Ifges(5,6)) * t97 + t173 * t95) * qJD(2) / 0.2e1 + t170 * t97 / 0.2e1) * t96 + ((t150 / 0.2e1 + t37 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1 + t22 * mrSges(5,2) - t30 * mrSges(4,3) - t15 * mrSges(6,3) + t120 * t96) * t97 + (t32 / 0.2e1 + t33 / 0.2e1 - t34 / 0.2e1 + t19 * mrSges(6,3) - t21 * mrSges(5,2) - t31 * mrSges(4,3) + (Ifges(6,6) / 0.2e1 + Ifges(4,6) / 0.2e1) * t98 - t119 * t96) * t95 + ((-t61 + t62) * t97 + (-t59 - t63) * t95 + m(5) * (-t21 * t95 + t22 * t97) + m(4) * (-t30 * t97 - t31 * t95)) * pkin(7)) * qJD(3); (t29 * t47 + t39 * t70 + t41 * t67) * t167 - 0.2e1 * pkin(2) * t51 + 0.2e1 * t65 * t49 + 0.2e1 * t29 * t69 + 0.2e1 * t47 * t50 + 0.2e1 * (m(5) * t65 + t68) * t40 + (t39 * t164 - t172 + t54) * t97 + (t41 * t164 + t170) * t95 + ((t164 * t67 + t74 + t75 + t76) * t97 + (0.2e1 * mrSges(6,3) * t70 + t71 + t72 - t73) * t95) * qJD(3); -t106 + (Ifges(6,3) * qJD(2) + (-t128 * t95 - t147 * t97) * qJD(3)) * t96 + (t26 + t28) * qJ(4) + (t58 + t63) * qJD(4) + m(6) * (qJ(4) * t2 + qJD(4) * t19 - t1 * t161) + m(5) * (-pkin(3) * t5 + qJ(4) * t4 + qJD(4) * t21) + (-Ifges(6,5) * t97 - t147 * t95) * t139 - t161 * t23 - pkin(3) * t25 - t7 * mrSges(4,2) + t8 * mrSges(4,1) - t1 * mrSges(6,1) + t2 * mrSges(6,2) + t4 * mrSges(5,3) - t5 * mrSges(5,1); m(6) * (qJ(4) * t39 + qJD(4) * t70 - t161 * t41) - t41 * mrSges(6,1) + t39 * mrSges(6,2) + t90 + t93 + t92 + qJD(4) * t107 + ((-pkin(3) * mrSges(5,2) + mrSges(6,3) * t161 - Ifges(6,5)) * t97 + (qJ(4) * t148 - t147) * t95 + (-m(5) * t109 + t118 + t68) * pkin(7)) * qJD(3); 0.2e1 * (mrSges(6,2) + mrSges(5,3) + (m(5) + m(6)) * qJ(4)) * qJD(4); -mrSges(5,2) * t126 + t78 + t79 + m(5) * t5 + m(6) * t1 + (-t129 + (-mrSges(5,1) - mrSges(6,1)) * t96) * qJD(2); m(6) * t41 + qJD(3) * t107; 0; 0; m(6) * t3 + t17; m(6) * t29 + t50; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t42(1), t42(2), t42(4), t42(7), t42(11); t42(2), t42(3), t42(5), t42(8), t42(12); t42(4), t42(5), t42(6), t42(9), t42(13); t42(7), t42(8), t42(9), t42(10), t42(14); t42(11), t42(12), t42(13), t42(14), t42(15);];
Mq = res;
