% Calculate time derivative of joint inertia matrix for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:14
% EndTime: 2019-03-08 19:51:17
% DurationCPUTime: 1.24s
% Computational Cost: add. (945->234), mult. (2295->346), div. (0->0), fcn. (1755->8), ass. (0->115)
t66 = cos(qJ(6));
t106 = qJD(6) * t66;
t67 = cos(qJ(4));
t109 = qJD(4) * t67;
t63 = sin(qJ(6));
t64 = sin(qJ(4));
t77 = t64 * t106 + t63 * t109;
t143 = 2 * qJ(3);
t142 = m(6) + m(7);
t43 = mrSges(7,1) * t63 + mrSges(7,2) * t66;
t140 = -mrSges(6,3) - t43;
t44 = -t64 * mrSges(6,2) - t67 * mrSges(6,3);
t45 = t64 * mrSges(5,1) + t67 * mrSges(5,2);
t139 = -t44 - t45;
t138 = -t64 * pkin(4) + qJ(5) * t67;
t61 = sin(pkin(6));
t68 = cos(qJ(2));
t127 = t61 * t68;
t62 = cos(pkin(6));
t23 = t127 * t67 + t62 * t64;
t65 = sin(qJ(2));
t128 = t61 * t65;
t99 = qJD(2) * t128;
t12 = qJD(4) * t23 - t64 * t99;
t11 = t64 * t12;
t102 = t64 * t127;
t13 = -qJD(4) * t102 + t109 * t62 - t67 * t99;
t129 = t13 * t67;
t24 = t62 * t67 - t102;
t75 = (t23 * t64 + t24 * t67) * qJD(4) - t11 - t129;
t14 = -t128 * t63 + t23 * t66;
t15 = t128 * t66 + t23 * t63;
t113 = qJD(2) * t68;
t98 = t61 * t113;
t3 = -qJD(6) * t15 + t13 * t66 - t63 * t98;
t4 = qJD(6) * t14 + t13 * t63 + t66 * t98;
t76 = -t66 * t3 - t63 * t4 + (t14 * t63 - t15 * t66) * qJD(6);
t137 = 2 * m(7);
t136 = -t63 / 0.2e1;
t135 = t63 / 0.2e1;
t134 = t66 / 0.2e1;
t70 = -pkin(2) - pkin(8);
t133 = pkin(5) - t70;
t132 = Ifges(7,4) * t63;
t131 = Ifges(7,4) * t66;
t130 = Ifges(7,6) * t67;
t5 = t24 * t12;
t85 = Ifges(7,1) * t63 + t131;
t21 = t67 * Ifges(7,5) + t64 * t85;
t126 = t63 * t21;
t47 = Ifges(7,1) * t66 - t132;
t125 = t63 * t47;
t69 = -pkin(4) - pkin(9);
t124 = t63 * t69;
t123 = t64 * mrSges(7,3);
t84 = Ifges(7,2) * t66 + t132;
t20 = t64 * t84 + t130;
t122 = t66 * t20;
t46 = -Ifges(7,2) * t63 + t131;
t121 = t66 * t46;
t120 = t66 * t69;
t119 = -mrSges(5,1) + mrSges(6,2);
t118 = t24 * t109 - t11;
t117 = qJ(3) * t98 + qJD(3) * t128;
t115 = t63 ^ 2 + t66 ^ 2;
t112 = qJD(4) * t63;
t111 = qJD(4) * t64;
t110 = qJD(4) * t66;
t108 = qJD(6) * t63;
t107 = qJD(6) * t64;
t105 = qJD(6) * t67;
t104 = qJ(5) * qJD(4);
t103 = 0.2e1 * t67;
t101 = mrSges(5,2) + t140;
t95 = t66 * t109;
t100 = t77 * Ifges(7,5) + Ifges(7,6) * t95;
t94 = pkin(4) * t109 + t64 * t104 + qJD(3);
t93 = m(7) * t115;
t92 = m(6) * t70 - mrSges(6,1);
t90 = qJD(4) * t133;
t16 = (pkin(9) * qJD(4) - qJD(5)) * t67 + t94;
t27 = t64 * t90;
t36 = qJ(3) - t138;
t25 = pkin(9) * t64 + t36;
t38 = t133 * t67;
t9 = -t25 * t63 + t38 * t66;
t1 = qJD(6) * t9 + t16 * t66 - t27 * t63;
t10 = t25 * t66 + t38 * t63;
t2 = -qJD(6) * t10 - t16 * t63 - t27 * t66;
t89 = t1 * t63 + t2 * t66;
t87 = t10 * t66 - t9 * t63;
t86 = mrSges(7,1) * t66 - mrSges(7,2) * t63;
t83 = -Ifges(7,5) * t63 - Ifges(7,6) * t66;
t79 = -qJ(5) * t12 + qJD(5) * t24;
t78 = -t107 * t63 + t95;
t17 = mrSges(7,2) * t111 + mrSges(7,3) * t78;
t18 = -mrSges(7,1) * t111 - mrSges(7,3) * t77;
t34 = mrSges(7,1) * t67 - t123 * t63;
t35 = -mrSges(7,2) * t67 + t123 * t66;
t73 = t106 * t35 - t108 * t34 + t63 * t17 + t66 * t18;
t72 = m(7) * t76;
t71 = m(6) * t75;
t37 = t133 * t64;
t33 = t85 * qJD(6);
t32 = t84 * qJD(6);
t31 = (mrSges(5,1) * t67 - mrSges(5,2) * t64) * qJD(4);
t30 = (-mrSges(6,2) * t67 + mrSges(6,3) * t64) * qJD(4);
t29 = t86 * qJD(6);
t28 = t67 * t90;
t26 = t86 * t64;
t22 = -qJD(5) * t67 + t94;
t8 = -mrSges(7,1) * t78 + mrSges(7,2) * t77;
t7 = t47 * t107 + (-Ifges(7,5) * t64 + t67 * t85) * qJD(4);
t6 = t46 * t107 + (-Ifges(7,6) * t64 + t84 * t67) * qJD(4);
t19 = [0.2e1 * m(7) * (t14 * t3 + t15 * t4 - t5) + 0.2e1 * (m(5) + m(6)) * (t61 ^ 2 * t65 * t113 + t13 * t23 - t5); t12 * t26 + t14 * t18 + t15 * t17 + t24 * t8 + t3 * t34 + t4 * t35 + ((t30 + t31) * t65 + ((-mrSges(3,1) + mrSges(4,2)) * t65 + (-mrSges(3,2) + mrSges(4,3) - t139) * t68) * qJD(2)) * t61 + m(6) * (t22 * t128 + t36 * t98) + m(4) * (-pkin(2) * t99 + t117) + m(7) * (t1 * t15 + t10 * t4 + t12 * t37 + t14 * t2 - t24 * t28 + t3 * t9) + t71 * t70 - (mrSges(5,3) + mrSges(6,1)) * t75 + (t75 * t70 + t117) * m(5); (t1 * t10 + t2 * t9 + t28 * t37) * t137 + 0.2e1 * t28 * t26 + t31 * t143 + 0.2e1 * t2 * t34 + 0.2e1 * t1 * t35 + 0.2e1 * t36 * t30 - 0.2e1 * t37 * t8 + 0.2e1 * t10 * t17 + 0.2e1 * t9 * t18 + ((t122 + t126 + (-Ifges(5,4) - Ifges(6,6)) * t103) * qJD(4) + t100) * t67 + (t66 * t6 + t63 * t7 + (t66 * t21 + (-t20 - t130) * t63) * qJD(6) + ((0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) + t83) * t64 + (Ifges(5,2) - Ifges(6,2) + Ifges(6,3) - Ifges(5,1) - Ifges(7,3)) * t103) * qJD(4)) * t64 + 0.2e1 * (m(6) * t36 + t44) * t22 + (0.2e1 * mrSges(4,3) + 0.2e1 * t45 + (m(4) + m(5)) * t143) * qJD(3); m(4) * t99 + m(7) * ((t14 * t66 + t15 * t63) * t111 + t76 * t67 + t118) + m(5) * (t111 * t23 + t118 - t129) + t71; (m(7) * (t10 * t112 + t110 * t9 - t28) + t35 * t112 + t34 * t110 + t8) * t64 + (m(7) * (-qJD(4) * t37 - t10 * t106 + t108 * t9 - t89) - qJD(4) * t26 - t73) * t67; (0.1e1 - t115) * t64 * t109 * t137; t24 * t29 + t119 * t13 + t101 * t12 + m(7) * t79 + m(6) * (-pkin(4) * t13 + t79) - t69 * t72 + t76 * mrSges(7,3); m(7) * (-qJ(5) * t28 - qJD(5) * t37 + t1 * t124 + t120 * t2) + t18 * t120 + t17 * t124 + t7 * t134 + t6 * t136 - t28 * t43 - qJD(5) * t26 - t37 * t29 + qJ(5) * t8 - t89 * mrSges(7,3) + (t92 * qJD(5) - t32 * t134 - t33 * t135) * t64 + (-t122 / 0.2e1 - t126 / 0.2e1 + t67 * t83 / 0.2e1 + (t134 * t47 + t136 * t46) * t64 - t87 * mrSges(7,3) + (m(7) * t87 - t63 * t34 + t66 * t35) * t69) * qJD(6) + ((pkin(4) * mrSges(6,1) - Ifges(7,5) * t66 / 0.2e1 + Ifges(7,6) * t135 - Ifges(5,5) + Ifges(6,4)) * t64 + (t121 / 0.2e1 + t125 / 0.2e1 - qJ(5) * mrSges(6,1) - Ifges(5,6) + Ifges(6,5)) * t67 + (m(6) * t138 + t139) * t70) * qJD(4); t64 * t29 + (-t101 * t67 + (-m(6) * pkin(4) - mrSges(7,3) * t115 + t69 * t93 + t119) * t64) * qJD(4) + t142 * (qJD(5) * t64 + t67 * t104); 0.2e1 * qJ(5) * t29 + t32 * t63 - t33 * t66 + (-t121 - t125) * qJD(6) + 0.2e1 * (t142 * qJ(5) - t140) * qJD(5); m(6) * t13 - t72; m(7) * (qJD(6) * t87 + t89) + t92 * t111 + t73; (m(6) + t93) * t111; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 + (-Ifges(7,6) * t108 - Ifges(7,3) * qJD(4)) * t64 + t100; (t105 * t66 - t111 * t63) * mrSges(7,2) + (t105 * t63 + t110 * t64) * mrSges(7,1); ((-mrSges(7,2) * t69 - Ifges(7,6)) * t66 + (-mrSges(7,1) * t69 - Ifges(7,5)) * t63) * qJD(6); -t43 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
