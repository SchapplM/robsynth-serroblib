% Calculate time derivative of joint inertia matrix for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:30
% EndTime: 2019-12-31 19:12:34
% DurationCPUTime: 1.54s
% Computational Cost: add. (1861->214), mult. (3824->333), div. (0->0), fcn. (3219->6), ass. (0->120)
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t60 = -mrSges(6,1) * t81 + mrSges(6,2) * t79;
t181 = t60 - mrSges(5,1);
t77 = t79 ^ 2;
t78 = t81 ^ 2;
t177 = t77 + t78;
t134 = qJD(5) * t79;
t160 = sin(qJ(3));
t161 = cos(qJ(4));
t103 = t161 * t160;
t162 = cos(qJ(3));
t111 = qJD(3) * t162;
t80 = sin(qJ(4));
t55 = t80 * t162 + t103;
t38 = -qJD(3) * t103 - t55 * qJD(4) - t80 * t111;
t137 = t81 * t38;
t104 = t161 * t162;
t54 = t80 * t160 - t104;
t93 = t54 * t134 + t137;
t82 = -pkin(1) - pkin(6);
t180 = -pkin(7) + t82;
t110 = qJD(3) * t160;
t39 = qJD(3) * t104 - t54 * qJD(4) - t80 * t110;
t179 = t39 * t55;
t147 = t54 * t38;
t178 = -Ifges(4,1) + Ifges(4,2);
t148 = t54 * mrSges(5,3);
t102 = mrSges(6,1) * t79 + mrSges(6,2) * t81;
t27 = t102 * t54;
t176 = -t27 - t148;
t175 = -mrSges(4,1) * t110 - mrSges(4,2) * t111;
t133 = qJD(5) * t81;
t141 = t79 * mrSges(6,3);
t32 = -mrSges(6,2) * t55 + t141 * t54;
t144 = t54 * t81;
t33 = mrSges(6,1) * t55 + mrSges(6,3) * t144;
t174 = -t33 * t133 - t32 * t134;
t153 = t38 * t79;
t94 = t54 * t133 - t153;
t70 = t160 * pkin(3) + qJ(2);
t28 = pkin(4) * t55 + pkin(8) * t54 + t70;
t59 = t180 * t160;
t90 = t180 * t162;
t43 = t161 * t59 + t80 * t90;
t13 = t28 * t81 - t43 * t79;
t14 = t28 * t79 + t43 * t81;
t173 = -t13 * t133 - t14 * t134;
t172 = t81 * t32 - t79 * t33;
t10 = -mrSges(6,2) * t39 + mrSges(6,3) * t94;
t9 = mrSges(6,1) * t39 - mrSges(6,3) * t93;
t171 = t81 * t10 - t79 * t9;
t170 = -t161 * t38 - t39 * t80;
t169 = 2 * qJD(2);
t168 = m(5) * pkin(3);
t56 = t102 * qJD(5);
t166 = pkin(4) * t56;
t63 = pkin(3) * t111 + qJD(2);
t12 = pkin(4) * t39 - pkin(8) * t38 + t63;
t135 = qJD(4) * t80;
t53 = t90 * qJD(3);
t86 = t161 * t90;
t88 = qJD(3) * t59;
t16 = qJD(4) * t86 - t135 * t59 + t161 * t53 - t80 * t88;
t2 = qJD(5) * t13 + t12 * t79 + t16 * t81;
t165 = t2 * t81;
t3 = -qJD(5) * t14 + t12 * t81 - t16 * t79;
t164 = t3 * t79;
t158 = Ifges(6,4) * t79;
t157 = Ifges(6,4) * t81;
t156 = Ifges(6,6) * t79;
t17 = qJD(4) * t43 + t161 * t88 + t80 * t53;
t42 = t80 * t59 - t86;
t155 = t17 * t42;
t154 = t38 * t42;
t152 = t39 * t79;
t150 = t39 * t81;
t149 = t42 * t80;
t146 = t54 * t79;
t145 = t54 * t80;
t143 = t55 * mrSges(5,3);
t128 = t161 * pkin(3);
t73 = -t128 - pkin(4);
t142 = t73 * t56;
t136 = pkin(3) * qJD(4);
t132 = -0.2e1 * t39 * mrSges(5,3);
t131 = m(6) * t136;
t130 = pkin(3) * t135;
t118 = t79 * t161;
t115 = t81 * t161;
t113 = t177 * t39;
t109 = mrSges(5,1) * t130;
t108 = t60 * t130;
t107 = qJD(4) * t128;
t101 = Ifges(6,1) * t81 - t158;
t100 = -Ifges(6,2) * t79 + t157;
t99 = t17 * t54 - t154;
t98 = mrSges(6,3) * t107;
t97 = mrSges(5,2) * t107;
t96 = t77 * t98;
t95 = t78 * t98;
t57 = t100 * qJD(5);
t58 = t101 * qJD(5);
t61 = Ifges(6,2) * t81 + t158;
t62 = Ifges(6,1) * t79 + t157;
t92 = t62 * t133 - t61 * t134 + t81 * t57 + t79 * t58;
t91 = -t181 * t38 + t54 * t56 + (mrSges(6,3) * t177 - mrSges(5,2)) * t39;
t89 = t177 * t161;
t87 = -t164 + (-t13 * t81 - t14 * t79) * qJD(5);
t85 = t93 * Ifges(6,5) + Ifges(6,6) * t94 + Ifges(6,3) * t39;
t20 = t55 * Ifges(6,6) - t100 * t54;
t21 = t55 * Ifges(6,5) - t101 * t54;
t6 = Ifges(6,4) * t93 + Ifges(6,2) * t94 + t39 * Ifges(6,6);
t7 = Ifges(6,1) * t93 + Ifges(6,4) * t94 + t39 * Ifges(6,5);
t74 = Ifges(6,5) * t133;
t84 = -t16 * mrSges(5,2) + mrSges(6,3) * t165 + t79 * t7 / 0.2e1 - t20 * t134 / 0.2e1 - t61 * t153 / 0.2e1 + t62 * t137 / 0.2e1 + t39 * (Ifges(6,5) * t79 + Ifges(6,6) * t81) / 0.2e1 + t42 * t56 - Ifges(5,6) * t39 + Ifges(5,5) * t38 - t58 * t144 / 0.2e1 + t55 * (-Ifges(6,6) * t134 + t74) / 0.2e1 + t81 * t6 / 0.2e1 + t181 * t17 + (qJD(5) * t62 + t57) * t146 / 0.2e1 + (t54 * t61 + t21) * t133 / 0.2e1;
t83 = m(6) * (-t164 + t165 + t173) + t171;
t72 = pkin(3) * t80 + pkin(8);
t8 = -mrSges(6,1) * t94 + mrSges(6,2) * t93;
t1 = [t43 * t132 + 0.2e1 * t63 * (mrSges(5,1) * t55 - mrSges(5,2) * t54) + 0.2e1 * t70 * (mrSges(5,1) * t39 + mrSges(5,2) * t38) + 0.2e1 * t42 * t8 + 0.2e1 * t2 * t32 + 0.2e1 * t3 * t33 + 0.2e1 * t13 * t9 + 0.2e1 * t14 * t10 - 0.2e1 * t16 * t143 - t7 * t144 + t55 * t85 + t6 * t146 + 0.2e1 * (-t38 * t55 + t39 * t54) * Ifges(5,4) + 0.2e1 * t176 * t17 + (mrSges(4,1) * t160 + mrSges(4,2) * t162 + mrSges(3,3)) * t169 + (0.2e1 * (mrSges(4,1) * t162 - mrSges(4,2) * t160) * qJD(3) + ((m(4) + m(3)) * t169)) * qJ(2) + 0.2e1 * Ifges(5,2) * t179 + (0.2e1 * Ifges(4,4) * t160 + t178 * t162) * t110 + (-0.2e1 * Ifges(4,4) * t162 + t178 * t160) * t111 - 0.2e1 * Ifges(5,1) * t147 + t93 * t21 + t94 * t20 + 0.2e1 * mrSges(5,3) * t154 + 0.2e1 * m(5) * (t16 * t43 + t63 * t70 + t155) + 0.2e1 * m(6) * (t13 * t3 + t14 * t2 + t155) + t39 * (Ifges(6,3) * t55 + (-Ifges(6,5) * t81 + t156) * t54); t54 * t8 + t172 * t39 + (t27 + 0.2e1 * t148) * t38 + m(6) * (-t13 * t152 + t14 * t150 + t99) + m(5) * (t39 * t43 + t99) + (t132 + (-t79 * t32 - t81 * t33) * qJD(5) + m(5) * t16 + t83) * t55; 0.2e1 * m(5) * (-t147 + t179) + 0.2e1 * m(6) * (t113 * t55 - t147); t73 * t8 - t3 * t141 + m(6) * (t73 * t17 + (t115 * t14 - t118 * t13 + t149) * t136) + t84 + (-t161 * t17 + t16 * t80 + (t161 * t43 + t149) * qJD(4)) * t168 - Ifges(4,5) * t110 - Ifges(4,6) * t111 + t175 * t82 + t170 * mrSges(5,3) * pkin(3) + t176 * t130 + t173 * mrSges(6,3) + (-t143 + t172) * t107 + (m(6) * (t87 + t165) + t171 + t174) * t72; ((t161 * t55 + t145) * qJD(4) - t170) * t168 + m(6) * (-t73 * t38 + t72 * t113 + (t55 * t89 + t145) * t136) + t91 + t175; 0.2e1 * t95 + 0.2e1 * t96 - 0.2e1 * t97 - 0.2e1 * t109 + 0.2e1 * (t72 * t89 + t73 * t80) * t131 + 0.2e1 * t108 + 0.2e1 * t142 + t92; (-m(6) * t17 - t8) * pkin(4) + (t83 + t174) * pkin(8) + t84 + t87 * mrSges(6,3); m(6) * (pkin(4) * t38 + pkin(8) * t113) + t91; t96 + t95 - t97 - t109 + (-pkin(4) * t80 + pkin(8) * t89) * t131 + t108 + t142 - t166 + t92; t92 - 0.2e1 * t166; mrSges(6,1) * t3 - mrSges(6,2) * t2 + t85; (t134 * t55 - t150) * mrSges(6,2) + (-t133 * t55 - t152) * mrSges(6,1); t74 + (-mrSges(6,1) * t118 - mrSges(6,2) * t115) * t136 + (t60 * t72 - t156) * qJD(5); t74 + (pkin(8) * t60 - t156) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
