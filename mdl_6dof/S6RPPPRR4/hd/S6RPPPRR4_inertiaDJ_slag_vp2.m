% Calculate time derivative of joint inertia matrix for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:18
% EndTime: 2019-03-09 01:35:20
% DurationCPUTime: 1.16s
% Computational Cost: add. (946->201), mult. (1869->317), div. (0->0), fcn. (1347->6), ass. (0->102)
t48 = sin(qJ(6));
t50 = cos(qJ(6));
t97 = t48 ^ 2 + t50 ^ 2;
t126 = m(7) * (-0.1e1 + t97);
t51 = cos(qJ(5));
t89 = qJD(6) * t51;
t49 = sin(qJ(5));
t94 = qJD(5) * t49;
t59 = t48 * t89 + t50 * t94;
t64 = mrSges(7,1) * t48 + mrSges(7,2) * t50;
t20 = t64 * t51;
t76 = t50 * t89;
t58 = t48 * t94 - t76;
t7 = mrSges(7,1) * t58 + mrSges(7,2) * t59;
t125 = t20 * t94 + t51 * t7;
t29 = -mrSges(7,1) * t50 + mrSges(7,2) * t48;
t122 = -mrSges(6,1) + t29;
t121 = Ifges(6,2) - Ifges(6,1);
t43 = t49 ^ 2;
t45 = t51 ^ 2;
t120 = t43 - t45;
t68 = m(7) * pkin(5) - t122;
t119 = 0.2e1 * m(7);
t117 = t48 / 0.2e1;
t116 = -t50 / 0.2e1;
t93 = qJD(5) * t51;
t78 = t49 * t93;
t115 = t78 * t126;
t114 = pkin(5) * t49;
t113 = pkin(5) * t51;
t112 = pkin(8) * t49;
t111 = pkin(8) * t51;
t107 = Ifges(7,4) * t48;
t106 = Ifges(7,4) * t50;
t105 = Ifges(7,6) * t48;
t104 = Ifges(7,6) * t49;
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t52 = -pkin(1) - pkin(2);
t98 = t47 * qJ(2) + t46 * t52;
t22 = -qJ(4) + t98;
t95 = qJD(2) * t47;
t36 = -qJD(4) + t95;
t103 = t22 * t36;
t102 = t48 * t49;
t101 = t49 * t50;
t99 = t51 * mrSges(7,3);
t27 = mrSges(7,2) * t49 + t48 * t99;
t100 = t50 * t27;
t92 = qJD(6) * t48;
t91 = qJD(6) * t49;
t90 = qJD(6) * t50;
t88 = t46 * qJD(2);
t87 = -0.2e1 * t88;
t85 = t43 * t88;
t84 = t45 * t95;
t83 = t49 * t88;
t71 = -t46 * qJ(2) + t47 * t52;
t67 = pkin(3) - t71;
t21 = pkin(7) + t67;
t82 = t21 * t93;
t80 = t48 * t93;
t35 = t45 * t88;
t73 = t50 * t93;
t13 = t111 + t22 - t114;
t57 = qJD(6) * t13 + t82 + t83;
t66 = -t21 * t91 + (-t112 - t113) * qJD(5) + t36;
t1 = t48 * t66 + t50 * t57;
t3 = -t102 * t21 + t13 * t50;
t72 = -qJD(6) * t3 + t1;
t70 = 0.2e1 * t47 * t93;
t65 = -t49 * mrSges(6,1) - t51 * mrSges(6,2);
t63 = Ifges(7,1) * t50 - t107;
t31 = Ifges(7,1) * t48 + t106;
t62 = -Ifges(7,2) * t48 + t106;
t30 = Ifges(7,2) * t50 + t107;
t18 = t102 * t47 + t46 * t50;
t60 = t101 * t47 - t46 * t48;
t61 = t18 * t48 + t50 * t60;
t56 = Ifges(7,6) * t76 + (-Ifges(7,6) * t102 - Ifges(7,3) * t51) * qJD(5) + t59 * Ifges(7,5);
t2 = -t48 * t57 + t50 * t66;
t4 = t101 * t21 + t13 * t48;
t55 = t1 * t50 - t2 * t48 - t3 * t90 - t4 * t92;
t8 = qJD(6) * t18 - t47 * t73;
t9 = qJD(6) * t60 + t47 * t80;
t54 = -t9 * t48 + t8 * t50 + (-t18 * t50 + t48 * t60) * qJD(6);
t11 = -mrSges(7,1) * t93 - mrSges(7,3) * t59;
t12 = mrSges(7,2) * t93 - mrSges(7,3) * t58;
t28 = -mrSges(7,1) * t49 + t50 * t99;
t53 = -qJD(5) * t20 - t48 * t11 + t50 * t12 - t27 * t92 - t28 * t90;
t39 = Ifges(7,5) * t90;
t38 = mrSges(6,2) * t94;
t26 = t63 * qJD(6);
t25 = t62 * qJD(6);
t24 = -mrSges(6,1) * t93 + t38;
t23 = t64 * qJD(6);
t16 = -Ifges(7,5) * t49 - t51 * t63;
t15 = -t51 * t62 - t104;
t14 = t21 * t35;
t6 = t31 * t89 + (-Ifges(7,5) * t51 + t49 * t63) * qJD(5);
t5 = t30 * t89 + (-Ifges(7,6) * t51 + t49 * t62) * qJD(5);
t10 = [mrSges(5,2) * t87 + 0.2e1 * t22 * t24 + 0.2e1 * t1 * t27 + 0.2e1 * t2 * t28 + 0.2e1 * t3 * t11 + 0.2e1 * t4 * t12 + 0.2e1 * m(5) * (t67 * t88 + t103) + 0.2e1 * m(6) * (t21 * t85 + t103 + t14) + 0.2e1 * mrSges(4,2) * t95 + (-t21 ^ 2 * t78 + t4 * t1 + t2 * t3 + t14) * t119 + 0.2e1 * (-mrSges(5,3) + t65) * t36 + t59 * t16 - t58 * t15 + (-t56 + (Ifges(7,3) + t121) * t93) * t49 + 0.2e1 * (m(4) * (-t46 * t71 + t47 * t98) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + (-t50 * t6 + t48 * t5 - (-Ifges(7,5) * t50 + t105) * t93 - t20 * t87 + t121 * t94) * t51 + 0.2e1 * (mrSges(4,1) + (t45 + t43) * mrSges(6,3)) * t88 - 0.2e1 * t125 * t21 + 0.2e1 * (t49 * t94 - t51 * t93) * Ifges(6,4); t18 * t11 - t60 * t12 + t8 * t27 + t9 * t28 + t125 * t47 + m(7) * (t21 * t49 * t70 - t1 * t60 + t18 * t2 + t3 * t9 + t4 * t8) + (t24 - m(7) * t84 + m(6) * (-t43 * t95 + t36 - t84) + m(5) * (t36 - t95)) * t46; (-t47 ^ 2 * t78 + t18 * t9 - t60 * t8) * t119; t49 * t7 + (m(7) * (-t4 * t101 + t3 * t102 + t120 * t21) - t49 * t100 + t28 * t102) * qJD(5) + (m(7) * (t55 - t83) + t53) * t51; m(7) * (t54 * t51 + (-t120 * t47 + t49 * t61) * qJD(5)); -0.2e1 * t115; m(5) * t88 + (-t7 + (-t48 * t28 + t100) * qJD(5)) * t51 + m(7) * (-t3 * t80 + t4 * t73 + t35) + m(6) * (t35 + t85) + (m(7) * (t55 - 0.2e1 * t82) + t53) * t49; m(7) * (-t61 * t93 + (t54 + t70) * t49); -t120 * qJD(5) * t126; 0.2e1 * t115; -pkin(5) * t7 + (-mrSges(6,2) * t88 - t39 / 0.2e1 + (-t21 * t68 + Ifges(6,5)) * qJD(5)) * t49 + (-t30 * t94 / 0.2e1 - t2 * mrSges(7,3) + t6 / 0.2e1 + (-t4 * mrSges(7,3) - t15 / 0.2e1 + t104 / 0.2e1) * qJD(6) + (m(7) * (-qJD(6) * t4 - t2) - qJD(6) * t27 - t11) * pkin(8)) * t48 + (t31 * t94 / 0.2e1 + qJD(6) * t16 / 0.2e1 + t5 / 0.2e1 + t72 * mrSges(7,3) + (m(7) * t72 - qJD(6) * t28 + t12) * pkin(8)) * t50 + (t26 * t116 + t25 * t117 - t21 * t23 + (t50 * t30 / 0.2e1 + t31 * t117) * qJD(6) + (-t21 * mrSges(6,2) - Ifges(7,5) * t48 / 0.2e1 + Ifges(7,6) * t116 + Ifges(6,6)) * qJD(5) + t68 * t88) * t51; (m(7) * pkin(8) + mrSges(7,3)) * t54 + ((qJD(5) * mrSges(6,2) + t23) * t51 + t68 * t94) * t47; t49 * t23 + t38 + (m(7) * (-t97 * t112 - t113) - t97 * t49 * mrSges(7,3) + t122 * t51) * qJD(5); -t51 * t23 + (t49 * t29 + m(7) * (t97 * t111 - t114) + t97 * t99 + t65) * qJD(5); -0.2e1 * pkin(5) * t23 + t25 * t50 + t26 * t48 + (-t30 * t48 + t31 * t50) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t56; mrSges(7,1) * t9 - mrSges(7,2) * t8; t7; (t48 * t91 - t73) * mrSges(7,2) + (-t49 * t90 - t80) * mrSges(7,1); t39 + (pkin(8) * t29 - t105) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
