% Calculate time derivative of joint inertia matrix for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:38
% EndTime: 2019-03-09 01:33:41
% DurationCPUTime: 1.25s
% Computational Cost: add. (1891->219), mult. (3618->337), div. (0->0), fcn. (3335->8), ass. (0->103)
t60 = sin(pkin(10));
t62 = cos(pkin(10));
t116 = cos(qJ(5));
t83 = qJD(5) * t116;
t65 = sin(qJ(5));
t96 = qJD(5) * t65;
t127 = -t60 * t96 + t62 * t83;
t103 = t65 * t62;
t39 = t116 * t60 + t103;
t66 = cos(qJ(6));
t64 = sin(qJ(6));
t94 = qJD(6) * t64;
t71 = -t127 * t66 + t39 * t94;
t61 = sin(pkin(9));
t63 = cos(pkin(9));
t67 = -pkin(1) - pkin(2);
t101 = t63 * qJ(2) + t61 * t67;
t40 = -qJ(4) + t101;
t117 = pkin(7) - t40;
t29 = t117 * t62;
t97 = qJD(2) * t63;
t49 = -qJD(4) + t97;
t87 = t65 * t117;
t8 = -t29 * t83 + t49 * t103 + (qJD(5) * t87 + t116 * t49) * t60;
t129 = -0.2e1 * t8;
t128 = t127 * t39;
t37 = t39 * qJD(5);
t20 = -t37 * mrSges(6,1) - mrSges(6,2) * t127;
t93 = qJD(6) * t66;
t72 = t127 * t64 + t39 * t93;
t126 = (t60 ^ 2 + t62 ^ 2) * t49;
t45 = -mrSges(7,1) * t66 + mrSges(7,2) * t64;
t125 = -m(7) * pkin(5) - mrSges(6,1) + t45;
t124 = 2 * mrSges(6,3);
t80 = t117 * t116;
t16 = -t65 * t29 - t60 * t80;
t123 = 0.2e1 * t16;
t121 = t39 / 0.2e1;
t114 = Ifges(7,4) * t66;
t47 = Ifges(7,1) * t64 + t114;
t120 = t47 / 0.2e1;
t119 = pkin(5) * t37;
t118 = t16 * t8;
t115 = Ifges(7,4) * t64;
t113 = Ifges(7,6) * t64;
t28 = t127 * t61;
t31 = t39 * t61;
t112 = t28 * t31;
t110 = t37 * Ifges(7,5);
t109 = t37 * Ifges(7,6);
t86 = t116 * t62;
t70 = -t65 * t60 + t86;
t108 = t70 * Ifges(7,6);
t107 = t70 * t37;
t106 = t39 * t64;
t105 = t39 * t66;
t17 = -t116 * t29 + t60 * t87;
t84 = -t61 * qJ(2) + t63 * t67;
t82 = pkin(3) - t84;
t35 = t62 * pkin(4) + t82;
t18 = pkin(5) * t70 + t39 * pkin(8) + t35;
t3 = -t17 * t64 + t18 * t66;
t100 = t3 * qJD(6);
t4 = t17 * t66 + t18 * t64;
t99 = t4 * qJD(6);
t98 = qJD(2) * t61;
t95 = qJD(6) * t39;
t85 = (-t64 ^ 2 - t66 ^ 2) * t127;
t79 = t16 * t28 + t31 * t8;
t78 = t16 * t37 - t70 * t8;
t77 = mrSges(7,1) * t64 + mrSges(7,2) * t66;
t76 = Ifges(7,1) * t66 - t115;
t75 = -Ifges(7,2) * t64 + t114;
t32 = t70 * t61;
t24 = -t32 * t64 - t63 * t66;
t25 = t32 * t66 - t63 * t64;
t74 = t24 * t64 - t25 * t66;
t73 = -t28 * t70 + t31 * t37;
t27 = t61 * t37;
t12 = qJD(6) * t24 - t27 * t66;
t13 = -qJD(6) * t25 + t27 * t64;
t69 = t12 * t66 - t13 * t64 + (-t24 * t66 - t25 * t64) * qJD(6);
t68 = t71 * Ifges(7,5) + Ifges(7,6) * t72 - Ifges(7,3) * t37;
t9 = -mrSges(7,1) * t72 + mrSges(7,2) * t71;
t52 = Ifges(7,5) * t93;
t46 = Ifges(7,2) * t66 + t115;
t43 = t76 * qJD(6);
t42 = t75 * qJD(6);
t41 = t77 * qJD(6);
t23 = mrSges(7,1) * t70 + mrSges(7,3) * t105;
t22 = -mrSges(7,2) * t70 + mrSges(7,3) * t106;
t21 = t77 * t39;
t19 = pkin(8) * t127 - t119 + t98;
t15 = Ifges(7,5) * t70 - t39 * t76;
t14 = -t39 * t75 + t108;
t11 = mrSges(7,2) * t37 + mrSges(7,3) * t72;
t10 = -mrSges(7,1) * t37 - mrSges(7,3) * t71;
t7 = t49 * t86 + t29 * t96 + (qJD(5) * t80 - t65 * t49) * t60;
t6 = Ifges(7,1) * t71 + Ifges(7,4) * t72 - t110;
t5 = Ifges(7,4) * t71 + Ifges(7,2) * t72 - t109;
t2 = t19 * t66 - t64 * t7 - t99;
t1 = t19 * t64 + t66 * t7 + t100;
t26 = [0.2e1 * t35 * t20 + t21 * t129 + 0.2e1 * t1 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t3 * t10 + 0.2e1 * t4 * t11 + t5 * t106 + t9 * t123 + 0.2e1 * m(7) * (t1 * t4 + t2 * t3 + t118) + 0.2e1 * m(6) * (t17 * t7 + t35 * t98 + t118) - t6 * t105 + 0.2e1 * mrSges(4,2) * t97 + 0.2e1 * m(5) * (t126 * t40 + t82 * t98) + 0.2e1 * Ifges(6,1) * t128 - 0.2e1 * mrSges(5,3) * t126 - (t7 * t124 - t68) * t70 + t71 * t15 + t72 * t14 + (-t123 * t127 + t129 * t39) * mrSges(6,3) + 0.2e1 * (mrSges(5,1) * t62 + mrSges(6,1) * t70 - mrSges(5,2) * t60 - mrSges(6,2) * t39 + mrSges(4,1)) * t98 + (t17 * t124 - (-Ifges(7,5) * t66 + t113) * t39 - ((2 * Ifges(6,2)) + Ifges(7,3)) * t70) * t37 + 0.2e1 * (m(4) * (t101 * t63 - t61 * t84) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (t127 * t70 - t37 * t39) * Ifges(6,4); t24 * t10 + t25 * t11 + t12 * t22 + t13 * t23 - t63 * t20 - t28 * t21 + t31 * t9 + m(7) * (t1 * t25 + t12 * t4 + t13 * t3 + t2 * t24 + t79) + m(5) * (-t97 + t126) * t61 + (-t127 * t31 + t27 * t70 - t28 * t39 + t32 * t37) * mrSges(6,3) + (-t17 * t27 + t32 * t7 - t61 * t97 + t79) * m(6); 0.2e1 * m(6) * (-t27 * t32 + t112) + 0.2e1 * m(7) * (t12 * t25 + t13 * t24 + t112); -t37 * t21 - t70 * t9 + m(7) * t78 + m(6) * (t127 * t17 + t39 * t7 + t78) + (m(7) * (t1 * t39 + t127 * t4 - t3 * t95) + t127 * t22 + t39 * t11 - t23 * t95) * t66 + (m(7) * (-t127 * t3 - t2 * t39 - t4 * t95) - t22 * t95 - t127 * t23 - t39 * t10) * t64; m(6) * (t127 * t32 - t27 * t39 + t73) + m(7) * (-t127 * t74 + t39 * t69 + t73); 0.2e1 * m(6) * (-t107 + t128) + 0.2e1 * m(7) * (-t39 * t85 - t107); m(7) * (t1 * t64 + t2 * t66 + (-t3 * t64 + t4 * t66) * qJD(6)) + t22 * t93 + t64 * t11 - t23 * t94 + t66 * t10 + (m(5) + m(6)) * t98 + t20; m(7) * (-qJD(6) * t74 + t12 * t64 + t13 * t66); 0; 0; t70 * t52 / 0.2e1 + t16 * t41 - Ifges(6,5) * t127 + Ifges(6,6) * t37 - t7 * mrSges(6,2) - pkin(5) * t9 + t125 * t8 + (t42 * t121 + t127 * t46 / 0.2e1 - t2 * mrSges(7,3) + t6 / 0.2e1 - t110 / 0.2e1 + (-t14 / 0.2e1 - t108 / 0.2e1 + t39 * t120 - t4 * mrSges(7,3)) * qJD(6) + (m(7) * (-t2 - t99) - qJD(6) * t22 - t10) * pkin(8)) * t64 + (-t39 * t43 / 0.2e1 - t127 * t120 + t1 * mrSges(7,3) - t109 / 0.2e1 + t5 / 0.2e1 + (t15 / 0.2e1 + t46 * t121 - t3 * mrSges(7,3)) * qJD(6) + (m(7) * (t1 - t100) + t11 - qJD(6) * t23) * pkin(8)) * t66; t27 * mrSges(6,2) + t31 * t41 + (m(7) * pkin(8) + mrSges(7,3)) * t69 + t125 * t28; t37 * t45 - t70 * t41 + m(7) * (-pkin(8) * t85 - t119) - mrSges(7,3) * t85 + t20; 0; -0.2e1 * pkin(5) * t41 + t42 * t66 + t43 * t64 + (-t64 * t46 + t66 * t47) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t68; mrSges(7,1) * t13 - mrSges(7,2) * t12; t9; -t41; t52 + (pkin(8) * t45 - t113) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t26(1) t26(2) t26(4) t26(7) t26(11) t26(16); t26(2) t26(3) t26(5) t26(8) t26(12) t26(17); t26(4) t26(5) t26(6) t26(9) t26(13) t26(18); t26(7) t26(8) t26(9) t26(10) t26(14) t26(19); t26(11) t26(12) t26(13) t26(14) t26(15) t26(20); t26(16) t26(17) t26(18) t26(19) t26(20) t26(21);];
Mq  = res;
