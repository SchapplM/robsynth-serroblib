% Calculate time derivative of joint inertia matrix for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:58
% EndTime: 2019-12-31 19:01:02
% DurationCPUTime: 1.53s
% Computational Cost: add. (1860->201), mult. (4050->320), div. (0->0), fcn. (3456->8), ass. (0->96)
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t109 = t75 ^ 2 + t78 ^ 2;
t139 = qJD(3) + qJD(4);
t76 = sin(qJ(4));
t77 = sin(qJ(3));
t79 = cos(qJ(4));
t80 = cos(qJ(3));
t58 = t76 * t77 - t79 * t80;
t42 = t139 * t58;
t144 = t109 * t42;
t140 = t109 * t79;
t105 = qJD(5) * t78;
t121 = t42 * t75;
t59 = t76 * t80 + t79 * t77;
t143 = t59 * t105 - t121;
t100 = -cos(pkin(9)) * pkin(1) - pkin(2);
t63 = -pkin(3) * t80 + t100;
t27 = t58 * pkin(4) - t59 * pkin(8) + t63;
t67 = sin(pkin(9)) * pkin(1) + pkin(6);
t129 = pkin(7) + t67;
t56 = t129 * t80;
t99 = t129 * t77;
t29 = t79 * t56 - t76 * t99;
t12 = t27 * t78 - t29 * t75;
t28 = t76 * t56 + t79 * t99;
t94 = qJD(3) * t129;
t53 = t77 * t94;
t92 = t80 * t94;
t15 = -t28 * qJD(4) - t79 * t53 - t76 * t92;
t43 = t139 * t59;
t132 = pkin(4) * t43;
t17 = qJD(3) * t77 * pkin(3) + pkin(8) * t42 + t132;
t2 = qJD(5) * t12 + t15 * t78 + t17 * t75;
t13 = t27 * t75 + t29 * t78;
t3 = -qJD(5) * t13 - t15 * t75 + t17 * t78;
t141 = t2 * t78 - t3 * t75;
t106 = qJD(5) * t75;
t103 = t59 * t106;
t112 = t78 * t42;
t85 = t103 + t112;
t138 = 2 * m(6);
t137 = 0.2e1 * pkin(3);
t16 = t29 * qJD(4) - t76 * t53 + t79 * t92;
t136 = 0.2e1 * t16;
t135 = 0.2e1 * t63;
t134 = m(5) / 0.2e1;
t128 = Ifges(6,4) * t75;
t127 = Ifges(6,4) * t78;
t126 = Ifges(6,6) * t75;
t125 = t16 * t28;
t120 = t58 * t43;
t119 = t58 * t76;
t118 = t59 * t75;
t117 = t59 * t78;
t114 = t76 * mrSges(5,1);
t64 = -mrSges(6,1) * t78 + mrSges(6,2) * t75;
t113 = t76 * t64;
t111 = t79 * mrSges(5,2);
t110 = -Ifges(6,5) * t112 + Ifges(6,3) * t43;
t108 = pkin(3) * qJD(4);
t107 = qJD(5) * t59;
t104 = 0.2e1 * t80;
t8 = -mrSges(6,1) * t143 + t85 * mrSges(6,2);
t101 = m(6) * t16 - t8;
t97 = -t106 / 0.2e1;
t96 = t43 * mrSges(5,1) - t42 * mrSges(5,2);
t95 = -(2 * Ifges(5,4)) - t126;
t93 = mrSges(6,3) * t140;
t91 = mrSges(6,1) * t75 + mrSges(6,2) * t78;
t90 = Ifges(6,1) * t78 - t128;
t89 = -Ifges(6,2) * t75 + t127;
t88 = Ifges(6,5) * t75 + Ifges(6,6) * t78;
t87 = t16 * t58 + t28 * t43;
t61 = t89 * qJD(5);
t62 = t90 * qJD(5);
t65 = Ifges(6,2) * t78 + t128;
t66 = Ifges(6,1) * t75 + t127;
t84 = t66 * t105 - t65 * t106 + t78 * t61 + t75 * t62;
t60 = t91 * qJD(5);
t83 = -mrSges(6,3) * t144 + t43 * t64 + t58 * t60 - t96;
t10 = mrSges(6,1) * t43 + t85 * mrSges(6,3);
t11 = -mrSges(6,2) * t43 - mrSges(6,3) * t143;
t36 = -mrSges(6,2) * t58 - mrSges(6,3) * t118;
t37 = mrSges(6,1) * t58 - mrSges(6,3) * t117;
t82 = m(6) * (-t12 * t105 - t13 * t106 + t141) + t78 * t11 - t75 * t10 - t37 * t105 - t36 * t106;
t20 = t58 * Ifges(6,6) + t89 * t59;
t21 = t58 * Ifges(6,5) + t90 * t59;
t6 = -t85 * Ifges(6,4) - Ifges(6,2) * t143 + t43 * Ifges(6,6);
t7 = -t85 * Ifges(6,1) - Ifges(6,4) * t143 + t43 * Ifges(6,5);
t71 = Ifges(6,5) * t105;
t81 = t20 * t97 + t21 * t105 / 0.2e1 + t28 * t60 + t75 * t7 / 0.2e1 - Ifges(5,5) * t42 - t61 * t118 / 0.2e1 + t62 * t117 / 0.2e1 + t58 * (-Ifges(6,6) * t106 + t71) / 0.2e1 + t78 * t6 / 0.2e1 - t15 * mrSges(5,2) + (-t112 / 0.2e1 + t59 * t97) * t66 + (t88 / 0.2e1 - Ifges(5,6)) * t43 + (t64 - mrSges(5,1)) * t16 - t143 * t65 / 0.2e1 + ((-t12 * t78 - t13 * t75) * qJD(5) + t141) * mrSges(6,3);
t70 = -pkin(3) * t79 - pkin(4);
t69 = pkin(3) * t76 + pkin(8);
t30 = t91 * t59;
t1 = [-t21 * t112 + t96 * t135 - 0.2e1 * t28 * t8 + t30 * t136 + 0.2e1 * t2 * t36 + 0.2e1 * t3 * t37 + 0.2e1 * t12 * t10 + 0.2e1 * t13 * t11 + t20 * t121 + 0.2e1 * (-t28 * t42 - t29 * t43) * mrSges(5,3) + (t12 * t3 + t13 * t2 + t125) * t138 + 0.2e1 * m(5) * (t15 * t29 + t125) + (-0.2e1 * t15 * mrSges(5,3) + ((2 * Ifges(5,2)) + Ifges(6,3)) * t43 - t95 * t42 + t110) * t58 + (mrSges(5,3) * t136 - 0.2e1 * Ifges(5,1) * t42 - t75 * t6 + t78 * t7 + (Ifges(6,5) * t78 + t95) * t43 + (-t78 * t20 - t75 * t21 - t58 * t88) * qJD(5)) * t59 + ((t100 * mrSges(4,2) + Ifges(4,4) * t80) * t104 + (m(5) * pkin(3) * t135 + 0.2e1 * t100 * mrSges(4,1) + (mrSges(5,1) * t58 + mrSges(5,2) * t59) * t137 - 0.2e1 * Ifges(4,4) * t77 + (-Ifges(4,2) + Ifges(4,1)) * t104) * t77) * qJD(3); t43 * t30 - t58 * t8 + m(6) * t87 + m(5) * (t15 * t59 - t29 * t42 + t87) + (-t42 * t36 + t59 * t11 - t37 * t107 + m(6) * (-t12 * t107 - t13 * t42 + t2 * t59)) * t78 + (-t36 * t107 + t42 * t37 - t59 * t10 + m(6) * (-t13 * t107 + t12 * t42 - t3 * t59)) * t75; 0.2e1 * m(6) * (-t144 * t59 + t120) + 0.2e1 * m(5) * (-t42 * t59 + t120); (Ifges(4,5) * t80 - Ifges(4,6) * t77 + (-mrSges(4,1) * t80 + mrSges(4,2) * t77) * t67) * qJD(3) + t101 * t70 + t82 * t69 + (m(5) * (t15 * t76 - t16 * t79) + (t79 * t42 - t76 * t43) * mrSges(5,3) + ((-t58 * mrSges(5,3) + t78 * t36 - t75 * t37 + m(5) * t29 + m(6) * (-t12 * t75 + t13 * t78)) * t79 + (t59 * mrSges(5,3) + t30 + (m(6) + m(5)) * t28) * t76) * qJD(4)) * pkin(3) + t81; m(6) * (-t144 * t69 + t43 * t70) + (-mrSges(4,1) * t77 - mrSges(4,2) * t80) * qJD(3) + ((-t42 * t76 - t43 * t79) * t134 + (m(6) * (t140 * t59 + t119) / 0.2e1 + (t59 * t79 + t119) * t134) * qJD(4)) * t137 + t83; 0.2e1 * t70 * t60 + (-0.2e1 * t111 - 0.2e1 * t114 + (t140 * t69 + t70 * t76) * t138 + 0.2e1 * t113 + 0.2e1 * t93) * t108 + t84; -t101 * pkin(4) + t82 * pkin(8) + t81; m(6) * (-pkin(8) * t144 - t132) + t83; (t70 - pkin(4)) * t60 + (-t111 - t114 + m(6) * (-pkin(4) * t76 + t140 * pkin(8)) + t113 + t93) * t108 + t84; -0.2e1 * pkin(4) * t60 + t84; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t103 - Ifges(6,6) * t143 + t110; t8; t71 - t91 * t79 * t108 + (t64 * t69 - t126) * qJD(5); t71 + (t64 * pkin(8) - t126) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
