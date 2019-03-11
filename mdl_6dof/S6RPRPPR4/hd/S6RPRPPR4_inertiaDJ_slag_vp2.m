% Calculate time derivative of joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:48
% EndTime: 2019-03-09 02:46:52
% DurationCPUTime: 2.00s
% Computational Cost: add. (2896->312), mult. (6671->468), div. (0->0), fcn. (6374->8), ass. (0->122)
t113 = cos(pkin(10));
t109 = t113 ^ 2;
t111 = sin(pkin(10));
t159 = qJD(4) * (t111 ^ 2 + t109);
t157 = Ifges(5,1) + Ifges(6,1);
t112 = sin(pkin(9));
t146 = pkin(7) + qJ(2);
t101 = t146 * t112;
t114 = cos(pkin(9));
t103 = t146 * t114;
t116 = sin(qJ(3));
t118 = cos(qJ(3));
t156 = -t118 * t101 - t103 * t116;
t155 = (-Ifges(5,6) + Ifges(6,6)) * t111 + (Ifges(6,4) + Ifges(5,5)) * t113;
t154 = 2 * m(5);
t153 = 2 * m(6);
t152 = 2 * m(7);
t151 = -2 * Ifges(4,4);
t149 = pkin(4) + pkin(5);
t73 = -t116 * t101 + t103 * t118;
t97 = t112 * t118 + t116 * t114;
t46 = qJD(2) * t97 + qJD(3) * t73;
t147 = t46 * t156;
t145 = -pkin(8) + qJ(4);
t95 = t112 * t116 - t118 * t114;
t88 = t95 * qJD(3);
t89 = t97 * qJD(3);
t40 = pkin(3) * t89 + qJ(4) * t88 - qJD(4) * t97;
t45 = -t95 * qJD(2) + qJD(3) * t156;
t18 = t111 * t40 + t113 * t45;
t137 = t111 * t88;
t52 = -mrSges(5,2) * t89 + mrSges(5,3) * t137;
t55 = mrSges(6,2) * t137 + mrSges(6,3) * t89;
t144 = t52 + t55;
t135 = t113 * t88;
t53 = mrSges(5,1) * t89 + mrSges(5,3) * t135;
t54 = -t89 * mrSges(6,1) - mrSges(6,2) * t135;
t143 = t53 - t54;
t127 = -pkin(2) * t114 - pkin(1);
t60 = pkin(3) * t95 - qJ(4) * t97 + t127;
t31 = t111 * t60 + t113 * t73;
t48 = -mrSges(5,1) * t137 - mrSges(5,2) * t135;
t115 = sin(qJ(6));
t117 = cos(qJ(6));
t120 = t111 * t115 + t113 * t117;
t86 = t120 * qJD(6);
t96 = t111 * t117 - t113 * t115;
t87 = t96 * qJD(6);
t142 = -Ifges(7,5) * t86 - Ifges(7,6) * t87;
t141 = Ifges(5,4) * t111;
t140 = Ifges(5,4) * t113;
t139 = Ifges(6,5) * t111;
t138 = Ifges(6,5) * t113;
t136 = t111 * t97;
t134 = t113 * t97;
t133 = qJ(5) * t113;
t131 = qJ(4) * t159;
t129 = qJD(5) * t111;
t26 = t95 * qJ(5) + t31;
t49 = t96 * t97;
t24 = qJD(6) * t49 - t120 * t88;
t25 = -t86 * t97 - t88 * t96;
t128 = -Ifges(7,5) * t24 - Ifges(7,6) * t25 + Ifges(7,3) * t89;
t126 = t89 * mrSges(4,1) - t88 * mrSges(4,2);
t57 = t87 * mrSges(7,1) - mrSges(7,2) * t86;
t123 = qJ(5) * t111 + pkin(3);
t41 = t111 * t45;
t17 = t113 * t40 - t41;
t65 = t111 * t73;
t30 = t113 * t60 - t65;
t10 = t89 * qJ(5) + t95 * qJD(5) + t18;
t47 = -mrSges(6,1) * t137 + mrSges(6,3) * t135;
t7 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t13 = t65 - t149 * t95 + (-pkin(8) * t97 - t60) * t113;
t19 = pkin(8) * t136 + t26;
t3 = -t115 * t19 + t117 * t13;
t4 = t115 * t13 + t117 * t19;
t102 = t145 * t113;
t99 = t145 * t111;
t72 = t102 * t117 + t115 * t99;
t70 = -t102 * t115 + t117 * t99;
t119 = -qJD(5) * t134 + t88 * t133 + t46;
t100 = -mrSges(6,1) * t113 - mrSges(6,3) * t111;
t98 = -pkin(4) * t113 - t123;
t90 = t113 * t149 + t123;
t69 = Ifges(7,1) * t96 - Ifges(7,4) * t120;
t68 = Ifges(7,4) * t96 - Ifges(7,2) * t120;
t67 = mrSges(7,1) * t120 + mrSges(7,2) * t96;
t64 = -mrSges(6,2) * t136 + mrSges(6,3) * t95;
t63 = -mrSges(6,1) * t95 + mrSges(6,2) * t134;
t62 = mrSges(5,1) * t95 - mrSges(5,3) * t134;
t61 = -mrSges(5,2) * t95 - mrSges(5,3) * t136;
t59 = -Ifges(7,1) * t86 - Ifges(7,4) * t87;
t58 = -Ifges(7,4) * t86 - Ifges(7,2) * t87;
t56 = (mrSges(6,1) * t111 - mrSges(6,3) * t113) * t97;
t50 = t120 * t97;
t44 = qJD(4) * t96 - qJD(6) * t72;
t43 = qJD(4) * t120 + qJD(6) * t70;
t38 = -mrSges(7,1) * t95 - mrSges(7,3) * t50;
t37 = mrSges(7,2) * t95 + mrSges(7,3) * t49;
t36 = t89 * Ifges(5,5) - (Ifges(5,1) * t113 - t141) * t88;
t35 = t89 * Ifges(6,4) - (Ifges(6,1) * t113 + t139) * t88;
t34 = t89 * Ifges(5,6) - (-Ifges(5,2) * t111 + t140) * t88;
t33 = t89 * Ifges(6,6) - (Ifges(6,3) * t111 + t138) * t88;
t32 = (pkin(4) * t111 - t133) * t97 - t156;
t29 = (-t111 * t149 + t133) * t97 + t156;
t28 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t27 = -pkin(4) * t95 - t30;
t21 = Ifges(7,1) * t50 + Ifges(7,4) * t49 - Ifges(7,5) * t95;
t20 = Ifges(7,4) * t50 + Ifges(7,2) * t49 - Ifges(7,6) * t95;
t16 = -pkin(4) * t137 + t119;
t15 = mrSges(7,2) * t89 + mrSges(7,3) * t25;
t14 = -mrSges(7,1) * t89 - mrSges(7,3) * t24;
t12 = -pkin(4) * t89 - t17;
t11 = -t137 * t149 + t119;
t9 = -pkin(8) * t137 + t10;
t8 = t41 - t149 * t89 + (pkin(8) * t88 - t40) * t113;
t6 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - t89 * Ifges(7,5);
t5 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - t89 * Ifges(7,6);
t2 = -qJD(6) * t4 - t115 * t9 + t117 * t8;
t1 = qJD(6) * t3 + t115 * t8 + t117 * t9;
t22 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t112 ^ 2 + t114 ^ 2) * qJD(2) - t89 * (Ifges(7,5) * t50 + Ifges(7,6) * t49) + 0.2e1 * (t156 * t88 - t73 * t89) * mrSges(4,3) - 0.2e1 * t156 * t48 + 0.2e1 * m(4) * (t45 * t73 - t147) + (t17 * t30 + t18 * t31 - t147) * t154 + (t10 * t26 + t12 * t27 + t16 * t32) * t153 + (t1 * t4 - t11 * t29 + t2 * t3) * t152 + (-0.2e1 * mrSges(4,3) * t45 + ((2 * Ifges(4,2)) + (2 * Ifges(6,2)) + (2 * Ifges(5,3)) + Ifges(7,3)) * t89 + t128 + (-t151 - 0.2e1 * t155) * t88) * t95 + 0.2e1 * t3 * t14 + 0.2e1 * t4 * t15 + t24 * t21 + t25 * t20 - 0.2e1 * t11 * t28 + 0.2e1 * t29 * t7 + 0.2e1 * t1 * t37 + 0.2e1 * t2 * t38 + 0.2e1 * t32 * t47 + t49 * t5 + t50 * t6 + 0.2e1 * t31 * t52 + 0.2e1 * t30 * t53 + 0.2e1 * t27 * t54 + 0.2e1 * t26 * t55 + 0.2e1 * t16 * t56 + 0.2e1 * t18 * t61 + 0.2e1 * t17 * t62 + 0.2e1 * t12 * t63 + 0.2e1 * t10 * t64 + 0.2e1 * t127 * t126 + (t111 * t33 - t111 * t34 + t113 * t35 + t113 * t36 + (t151 + t155) * t89 - ((2 * Ifges(4,1)) + t157 * t109 + ((Ifges(6,3) + Ifges(5,2)) * t111 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t113) * t111) * t88 + 0.2e1 * (mrSges(5,1) * t111 + mrSges(5,2) * t113 + mrSges(4,3)) * t46) * t97; -t120 * t14 + t96 * t15 - t86 * t37 - t87 * t38 + t143 * t113 + t144 * t111 + m(7) * (t1 * t96 - t120 * t2 - t3 * t87 - t4 * t86) + m(6) * (t10 * t111 - t113 * t12) + m(5) * (t111 * t18 + t113 * t17) + t126; (t120 * t87 - t86 * t96) * t152; -t95 * t142 / 0.2e1 + (-t17 * mrSges(5,3) + t12 * mrSges(6,2) + t35 / 0.2e1 + t36 / 0.2e1 + t46 * mrSges(5,2) + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t89 + (-t56 + t28) * qJD(5) + (-t62 + t63) * qJD(4) - t143 * qJ(4)) * t111 + (t18 * mrSges(5,3) + t10 * mrSges(6,2) - t33 / 0.2e1 + t34 / 0.2e1 - t46 * mrSges(5,1) + (Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t89 + (t61 + t64) * qJD(4) + t144 * qJ(4)) * t113 + m(7) * (t1 * t72 - t11 * t90 + t129 * t29 + t2 * t70 + t3 * t44 + t4 * t43) + m(5) * (-pkin(3) * t46 + (-t111 * t30 + t113 * t31) * qJD(4) + (-t17 * t111 + t18 * t113) * qJ(4)) + (-t1 * t120 - t2 * t96 + t3 * t86 - t4 * t87) * mrSges(7,3) - t89 * (Ifges(7,5) * t96 - Ifges(7,6) * t120) / 0.2e1 - t120 * t5 / 0.2e1 + m(6) * (t16 * t98 + (qJ(4) * t10 + qJD(4) * t26) * t113 + (qJ(4) * t12 + qJD(4) * t27 - qJD(5) * t32) * t111) + t43 * t37 + t44 * t38 - t45 * mrSges(4,2) - t46 * mrSges(4,1) - pkin(3) * t48 + t29 * t57 + t49 * t58 / 0.2e1 + t50 * t59 / 0.2e1 - t11 * t67 + t25 * t68 / 0.2e1 + t24 * t69 / 0.2e1 + t70 * t14 + t72 * t15 - t86 * t21 / 0.2e1 - t87 * t20 / 0.2e1 - Ifges(4,6) * t89 + t90 * t7 + t96 * t6 / 0.2e1 + t98 * t47 + t16 * t100 - (-t111 * (Ifges(5,2) * t113 + t141) / 0.2e1 + t111 * (-Ifges(6,3) * t113 + t139) / 0.2e1 + Ifges(4,5) + (t157 * t111 - t138 + t140) * t113 / 0.2e1) * t88; m(7) * (-t120 * t44 + t43 * t96 - t70 * t87 - t72 * t86); 0.2e1 * t90 * t57 - t120 * t58 + t96 * t59 - t87 * t68 - t86 * t69 + (t129 * t90 + t43 * t72 + t44 * t70) * t152 + t131 * t154 + (-t98 * t129 + t131) * t153 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * t159 + 0.2e1 * (-t100 + t67) * t129 + 0.2e1 * (-t120 * t43 - t44 * t96 + t70 * t86 - t72 * t87) * mrSges(7,3); m(5) * t46 + m(6) * t16 + m(7) * t11 + t47 + t48 - t7; 0; (-m(6) - m(7)) * t129 - t57; 0; t115 * t15 + t117 * t14 + (-t115 * t38 + t117 * t37) * qJD(6) + m(7) * (t1 * t115 + t117 * t2 + (-t115 * t3 + t117 * t4) * qJD(6)) + m(6) * t12 + t54; m(7) * (-t115 * t86 - t117 * t87 + (t115 * t120 + t117 * t96) * qJD(6)); m(7) * (t115 * t43 + t117 * t44 + (-t115 * t70 + t117 * t72) * qJD(6)) + m(6) * t111 * qJD(4) + (-t115 * t87 + t117 * t86 + (t115 * t96 - t117 * t120) * qJD(6)) * mrSges(7,3); 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - t128; -t57; mrSges(7,1) * t44 - mrSges(7,2) * t43 + t142; 0; (-mrSges(7,1) * t115 - mrSges(7,2) * t117) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
