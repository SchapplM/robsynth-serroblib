% Calculate time derivative of joint inertia matrix for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:55
% EndTime: 2019-12-31 20:16:59
% DurationCPUTime: 1.49s
% Computational Cost: add. (3423->218), mult. (7467->340), div. (0->0), fcn. (7356->8), ass. (0->108)
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t81 = -mrSges(6,1) * t99 + mrSges(6,2) * t96;
t155 = -mrSges(5,1) + t81;
t123 = qJD(5) * t99;
t141 = cos(qJ(4));
t100 = cos(qJ(2));
t95 = cos(pkin(9));
t125 = t100 * t95;
t94 = sin(pkin(9));
t98 = sin(qJ(2));
t73 = -t94 * t98 + t125;
t74 = t100 * t94 + t95 * t98;
t97 = sin(qJ(4));
t105 = t141 * t73 - t97 * t74;
t68 = t74 * qJD(2);
t69 = t73 * qJD(2);
t35 = t105 * qJD(4) + t141 * t69 - t97 * t68;
t52 = t141 * t74 + t97 * t73;
t107 = t52 * t123 + t35 * t96;
t128 = -qJ(3) - pkin(6);
t82 = t128 * t98;
t75 = t94 * t82;
t47 = -t74 * qJD(3) + (t128 * t125 - t75) * qJD(2);
t103 = -t69 * pkin(7) + t47;
t83 = t128 * t100;
t53 = t95 * t82 + t83 * t94;
t49 = -pkin(7) * t74 + t53;
t54 = -t95 * t83 + t75;
t50 = pkin(7) * t73 + t54;
t154 = t141 * t49 - t50 * t97;
t113 = qJD(2) * t128;
t48 = t95 * (qJD(3) * t100 + t98 * t113) + t94 * (-t98 * qJD(3) + t100 * t113);
t40 = -pkin(7) * t68 + t48;
t11 = t154 * qJD(4) + t97 * t103 + t141 * t40;
t89 = -pkin(2) * t100 - pkin(1);
t56 = -t73 * pkin(3) + t89;
t23 = -pkin(4) * t105 - t52 * pkin(8) + t56;
t25 = t141 * t50 + t97 * t49;
t15 = t23 * t99 - t25 * t96;
t36 = t52 * qJD(4) + t141 * t68 + t97 * t69;
t91 = qJD(2) * t98 * pkin(2);
t55 = pkin(3) * t68 + t91;
t17 = pkin(4) * t36 - pkin(8) * t35 + t55;
t2 = t15 * qJD(5) + t11 * t99 + t17 * t96;
t16 = t23 * t96 + t25 * t99;
t3 = -t16 * qJD(5) - t11 * t96 + t17 * t99;
t157 = t2 * t99 - t3 * t96;
t144 = pkin(2) * t94;
t88 = pkin(2) * t95 + pkin(3);
t65 = t141 * t88 - t97 * t144;
t59 = t65 * qJD(4);
t156 = t59 * mrSges(5,2);
t66 = t141 * t144 + t97 * t88;
t118 = (t96 ^ 2 + t99 ^ 2) * t59;
t153 = -t15 * t96 + t16 * t99;
t152 = 2 * m(5);
t151 = 2 * m(6);
t150 = -2 * mrSges(5,3);
t12 = t25 * qJD(4) - t141 * t103 + t97 * t40;
t149 = 0.2e1 * t12;
t148 = -0.2e1 * t154;
t147 = 0.2e1 * t55;
t146 = 0.2e1 * t89;
t140 = Ifges(6,4) * t96;
t139 = Ifges(6,4) * t99;
t138 = Ifges(6,6) * t96;
t137 = t12 * t154;
t60 = t66 * qJD(4);
t134 = t154 * t60;
t132 = t35 * t99;
t130 = t52 * t96;
t129 = t52 * t99;
t127 = Ifges(6,5) * t132 + Ifges(6,3) * t36;
t124 = qJD(5) * t96;
t122 = 0.2e1 * t100;
t121 = t52 * t124;
t116 = t68 * mrSges(4,1) + t69 * mrSges(4,2);
t115 = t36 * mrSges(5,1) + t35 * mrSges(5,2);
t114 = -(2 * Ifges(5,4)) - t138;
t112 = mrSges(6,1) * t96 + mrSges(6,2) * t99;
t111 = Ifges(6,1) * t99 - t140;
t110 = -Ifges(6,2) * t96 + t139;
t109 = Ifges(6,5) * t96 + Ifges(6,6) * t99;
t37 = mrSges(6,2) * t105 - mrSges(6,3) * t130;
t38 = -mrSges(6,1) * t105 - mrSges(6,3) * t129;
t108 = t99 * t37 - t96 * t38;
t106 = t121 - t132;
t79 = t110 * qJD(5);
t80 = t111 * qJD(5);
t84 = Ifges(6,2) * t99 + t140;
t85 = Ifges(6,1) * t96 + t139;
t104 = t85 * t123 - t84 * t124 + t99 * t79 + t96 * t80;
t13 = mrSges(6,1) * t36 + t106 * mrSges(6,3);
t14 = -mrSges(6,2) * t36 - t107 * mrSges(6,3);
t102 = -t96 * t13 + m(6) * (-t15 * t123 - t16 * t124 + t157) + t99 * t14 - t38 * t123 - t37 * t124;
t20 = -Ifges(6,6) * t105 + t110 * t52;
t21 = -Ifges(6,5) * t105 + t111 * t52;
t6 = -t106 * Ifges(6,4) - t107 * Ifges(6,2) + t36 * Ifges(6,6);
t7 = -t106 * Ifges(6,1) - t107 * Ifges(6,4) + t36 * Ifges(6,5);
t78 = t112 * qJD(5);
t90 = Ifges(6,5) * t123;
t101 = t21 * t123 / 0.2e1 - t154 * t78 + t85 * t132 / 0.2e1 + Ifges(5,5) * t35 + t96 * t7 / 0.2e1 - t79 * t130 / 0.2e1 + t80 * t129 / 0.2e1 - t105 * (-Ifges(6,6) * t124 + t90) / 0.2e1 + t99 * t6 / 0.2e1 - t11 * mrSges(5,2) + (t109 / 0.2e1 - Ifges(5,6)) * t36 - t107 * t84 / 0.2e1 + t155 * t12 - (t52 * t85 + t20) * t124 / 0.2e1 + ((-t15 * t99 - t16 * t96) * qJD(5) + t157) * mrSges(6,3);
t63 = pkin(8) + t66;
t62 = -pkin(4) - t65;
t30 = t112 * t52;
t8 = t107 * mrSges(6,1) - t106 * mrSges(6,2);
t1 = [t25 * t36 * t150 + 0.2e1 * t15 * t13 + 0.2e1 * t16 * t14 + t8 * t148 + t30 * t149 + 0.2e1 * t2 * t37 + 0.2e1 * t3 * t38 + 0.2e1 * t56 * t115 - 0.2e1 * t73 * Ifges(4,2) * t68 + 0.2e1 * t74 * t69 * Ifges(4,1) + t116 * t146 + (mrSges(5,3) * t148 - t96 * t20 + t99 * t21) * t35 + 0.2e1 * m(4) * (t47 * t53 + t48 * t54) + (t11 * t25 + t55 * t56 - t137) * t152 + (t15 * t3 + t16 * t2 - t137) * t151 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t100) * t122 + (m(4) * pkin(2) * t146 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t73 + mrSges(4,2) * t74) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t98 + (Ifges(3,1) - Ifges(3,2)) * t122) * t98) * qJD(2) - (mrSges(5,1) * t147 + t11 * t150 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t36 + t114 * t35 + t127) * t105 + (mrSges(5,2) * t147 + mrSges(5,3) * t149 + 0.2e1 * Ifges(5,1) * t35 - t96 * t6 + t99 * t7 + (Ifges(6,5) * t99 + t114) * t36 + (t105 * t109 - t99 * t20 - t96 * t21) * qJD(5)) * t52 + 0.2e1 * (-t74 * t68 + t73 * t69) * Ifges(4,4) + 0.2e1 * (-t47 * t74 + t48 * t73 - t53 * t69 - t54 * t68) * mrSges(4,3); (Ifges(3,5) * t100 - Ifges(3,6) * t98 + (-mrSges(3,1) * t100 + mrSges(3,2) * t98) * pkin(6)) * qJD(2) + (m(4) * (t47 * t95 + t48 * t94) + (-t94 * t68 - t95 * t69) * mrSges(4,3)) * pkin(2) + t102 * t63 + m(6) * (t12 * t62 - t134) + m(5) * (t11 * t66 - t12 * t65 - t134) + (-t65 * t35 - t66 * t36 + t60 * t52) * mrSges(5,3) + t101 + t47 * mrSges(4,1) - t48 * mrSges(4,2) + t60 * t30 + t62 * t8 - Ifges(4,6) * t68 + Ifges(4,5) * t69 + (m(5) * t25 + m(6) * t153 + mrSges(5,3) * t105 + t108) * t59; 0.2e1 * t62 * t78 - 0.2e1 * t156 + (t63 * t118 + t60 * t62) * t151 + (t59 * t66 - t60 * t65) * t152 + t104 + 0.2e1 * t155 * t60 + 0.2e1 * mrSges(6,3) * t118; m(4) * t91 + t99 * t13 + t96 * t14 + t108 * qJD(5) + m(6) * (t153 * qJD(5) + t2 * t96 + t3 * t99) + m(5) * t55 + t115 + t116; 0; 0; (-m(6) * t12 - t8) * pkin(4) + t102 * pkin(8) + t101; -t156 + (-pkin(4) + t62) * t78 + t104 + (m(6) * pkin(8) + mrSges(6,3)) * t118 + (-m(6) * pkin(4) + t155) * t60; 0; -0.2e1 * pkin(4) * t78 + t104; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t121 - t107 * Ifges(6,6) + t127; t90 - t112 * t59 + (t81 * t63 - t138) * qJD(5); -t78; t90 + (t81 * pkin(8) - t138) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
