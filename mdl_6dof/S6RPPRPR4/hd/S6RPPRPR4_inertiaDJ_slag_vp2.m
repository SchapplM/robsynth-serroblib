% Calculate time derivative of joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:17
% EndTime: 2019-03-09 01:46:21
% DurationCPUTime: 2.03s
% Computational Cost: add. (1996->242), mult. (3903->381), div. (0->0), fcn. (3480->8), ass. (0->114)
t114 = cos(pkin(10));
t72 = cos(qJ(4));
t113 = sin(pkin(10));
t70 = sin(qJ(4));
t91 = t113 * t70;
t40 = -t114 * t72 + t91;
t69 = sin(qJ(6));
t108 = qJD(6) * t69;
t39 = t40 * qJD(4);
t93 = t114 * t70;
t41 = t113 * t72 + t93;
t71 = cos(qJ(6));
t77 = t41 * t108 + t71 * t39;
t149 = t39 * t41;
t148 = Ifges(5,1) - Ifges(5,2);
t116 = t70 ^ 2 + t72 ^ 2;
t38 = t41 * qJD(4);
t22 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t107 = qJD(6) * t71;
t78 = t41 * t107 - t69 * t39;
t110 = qJD(4) * t72;
t111 = qJD(4) * t70;
t147 = -mrSges(5,1) * t110 + mrSges(5,2) * t111;
t146 = -mrSges(5,1) * t70 - mrSges(5,2) * t72;
t67 = sin(pkin(9));
t59 = t67 * qJD(2);
t48 = -pkin(4) * t111 + t59;
t145 = -m(6) * t48 - t22;
t98 = t39 * (t69 ^ 2 + t71 ^ 2);
t141 = m(6) * pkin(4);
t144 = t113 * t141 - mrSges(6,2);
t50 = -mrSges(7,1) * t71 + mrSges(7,2) * t69;
t97 = t114 * pkin(4);
t56 = -t97 - pkin(5);
t143 = m(7) * t56 - t114 * t141 - mrSges(6,1) + t50;
t18 = -pkin(5) * t38 - pkin(8) * t39 + t48;
t68 = cos(pkin(9));
t73 = -pkin(1) - pkin(2);
t117 = t68 * qJ(2) + t67 * t73;
t43 = -pkin(7) + t117;
t115 = qJ(5) - t43;
t31 = t115 * t72;
t17 = -t114 * t31 + t115 * t91;
t94 = -t67 * qJ(2) + t68 * t73;
t42 = pkin(3) - t94;
t37 = t72 * pkin(4) + t42;
t19 = -t40 * pkin(5) + t41 * pkin(8) + t37;
t3 = -t17 * t69 + t19 * t71;
t112 = qJD(2) * t68;
t88 = -qJD(5) + t112;
t89 = qJD(4) * t115;
t21 = t70 * t89 + t72 * t88;
t76 = -t70 * t88 + t72 * t89;
t8 = t113 * t76 + t114 * t21;
t1 = qJD(6) * t3 + t18 * t69 + t71 * t8;
t140 = t1 * t71;
t16 = -t113 * t31 - t115 * t93;
t7 = t113 * t21 - t114 * t76;
t139 = t16 * t7;
t134 = Ifges(7,4) * t69;
t133 = Ifges(7,4) * t71;
t132 = Ifges(7,6) * t69;
t34 = t40 * t67;
t26 = t34 * t69 - t68 * t71;
t30 = t67 * t38;
t12 = qJD(6) * t26 - t30 * t71;
t131 = t12 * t71;
t29 = t39 * t67;
t33 = t41 * t67;
t130 = t29 * t33;
t128 = t38 * mrSges(6,3);
t127 = t38 * t40;
t126 = t39 * mrSges(6,3);
t124 = t41 * t69;
t123 = t41 * t71;
t122 = t69 * mrSges(7,3);
t51 = Ifges(7,2) * t71 + t134;
t120 = t69 * t51;
t52 = Ifges(7,1) * t69 + t133;
t118 = t71 * t52;
t109 = qJD(6) * t41;
t103 = mrSges(7,3) * t107;
t102 = mrSges(7,3) * t108;
t96 = t113 * pkin(4);
t85 = -t16 * t29 + t33 * t7;
t84 = t16 * t38 + t40 * t7;
t83 = mrSges(7,1) * t69 + mrSges(7,2) * t71;
t82 = Ifges(7,1) * t71 - t134;
t81 = -Ifges(7,2) * t69 + t133;
t4 = t17 * t71 + t19 * t69;
t27 = -t34 * t71 - t68 * t69;
t80 = -t26 * t69 + t27 * t71;
t79 = -t29 * t40 + t33 * t38;
t13 = -qJD(6) * t27 + t30 * t69;
t75 = t131 - t13 * t69 + (-t26 * t71 - t27 * t69) * qJD(6);
t74 = Ifges(7,5) * t77 + Ifges(7,6) * t78 - Ifges(7,3) * t38;
t9 = -mrSges(7,1) * t78 + mrSges(7,2) * t77;
t58 = Ifges(7,5) * t107;
t55 = t96 + pkin(8);
t47 = t82 * qJD(6);
t46 = t81 * qJD(6);
t45 = t146 * qJD(4);
t44 = t83 * qJD(6);
t25 = -mrSges(7,1) * t40 + mrSges(7,3) * t123;
t24 = mrSges(7,2) * t40 + t122 * t41;
t23 = t83 * t41;
t15 = -t40 * Ifges(7,5) - t41 * t82;
t14 = -t40 * Ifges(7,6) - t41 * t81;
t11 = mrSges(7,2) * t38 + mrSges(7,3) * t78;
t10 = -mrSges(7,1) * t38 - mrSges(7,3) * t77;
t6 = Ifges(7,1) * t77 + Ifges(7,4) * t78 - t38 * Ifges(7,5);
t5 = Ifges(7,4) * t77 + Ifges(7,2) * t78 - t38 * Ifges(7,6);
t2 = -qJD(6) * t4 + t18 * t71 - t69 * t8;
t20 = [-0.2e1 * Ifges(6,1) * t149 + (0.2e1 * Ifges(5,4) * t72 + t148 * t70) * t110 + (-0.2e1 * Ifges(5,4) * t70 + t148 * t72) * t111 + t77 * t15 + t78 * t14 + 0.2e1 * (mrSges(5,1) * t72 - mrSges(5,2) * t70 + mrSges(4,1)) * t59 + 0.2e1 * (m(5) * (t116 * t43 * t68 + t42 * t67) + m(4) * (t117 * t68 - t67 * t94) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + (-(-Ifges(7,5) * t71 + t132) * t41 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t40) * t38 + 0.2e1 * (t8 * t40 - t7 * t41) * mrSges(6,3) + 0.2e1 * m(6) * (t17 * t8 + t37 * t48 + t139) + 0.2e1 * m(7) * (t1 * t4 + t2 * t3 + t139) + 0.2e1 * (-t116 * mrSges(5,3) + mrSges(4,2)) * t112 + 0.2e1 * t17 * t128 - t6 * t123 - t40 * t74 + 0.2e1 * t42 * t45 + 0.2e1 * t48 * (-mrSges(6,1) * t40 - mrSges(6,2) * t41) + 0.2e1 * t37 * t22 - 0.2e1 * t7 * t23 + 0.2e1 * t1 * t24 + 0.2e1 * t2 * t25 + 0.2e1 * t3 * t10 + 0.2e1 * t4 * t11 + 0.2e1 * (-t38 * t41 + t39 * t40) * Ifges(6,4) + 0.2e1 * (t126 + t9) * t16 + t5 * t124; t26 * t10 + t27 * t11 + t12 * t24 + t13 * t25 + t29 * t23 + t33 * t9 + m(7) * (t1 * t27 + t12 * t4 + t13 * t3 + t2 * t26 + t85) + m(6) * (-t17 * t30 - t34 * t8 + t85) + (t29 * t41 - t30 * t40 + t33 * t39 - t34 * t38) * mrSges(6,3) + (-t45 + m(5) * (-0.1e1 + t116) * t59 + t145) * t68; 0.2e1 * m(7) * (t12 * t27 + t13 * t26 - t130) + 0.2e1 * m(6) * (t30 * t34 - t130); -t38 * t23 + t40 * t9 + m(7) * t84 + m(6) * (-t17 * t39 + t41 * t8 + t84) + (m(7) * (t1 * t41 - t109 * t3 - t39 * t4) - t39 * t24 + t41 * t11 - t25 * t109) * t71 + (m(7) * (-t109 * t4 - t2 * t41 + t3 * t39) - t24 * t109 + t39 * t25 - t41 * t10) * t69; m(7) * (-t39 * t80 + t41 * t75 + t79) + m(6) * (-t30 * t41 + t34 * t39 + t79); 0.2e1 * m(6) * (t127 - t149) + 0.2e1 * m(7) * (-t41 * t98 + t127); t146 * t112 + t147 * t43 + t144 * t8 + t143 * t7 + (t41 * t51 + t15) * t107 / 0.2e1 + (qJD(6) * t52 + t46) * t124 / 0.2e1 + (m(7) * (t140 - t2 * t69 + (-t3 * t71 - t4 * t69) * qJD(6)) - t25 * t107 - t24 * t108 + t71 * t11 - t69 * t10) * t55 - (-t118 / 0.2e1 + t120 / 0.2e1 - Ifges(6,5)) * t39 - t47 * t123 / 0.2e1 - t2 * t122 - Ifges(5,5) * t110 - t14 * t108 / 0.2e1 - t40 * (-Ifges(7,6) * t108 + t58) / 0.2e1 - t3 * t103 - t4 * t102 + t69 * t6 / 0.2e1 - t38 * (Ifges(7,5) * t69 + Ifges(7,6) * t71) / 0.2e1 + t71 * t5 / 0.2e1 + t56 * t9 + t16 * t44 + Ifges(6,6) * t38 - t97 * t126 + t96 * t128 + mrSges(7,3) * t140 + Ifges(5,6) * t111; m(7) * t75 * t55 + mrSges(7,3) * t131 - t27 * t102 - t26 * t103 - t13 * t122 - t143 * t29 - t144 * t30 + t147 * t67 + t33 * t44; -mrSges(5,2) * t110 - mrSges(5,1) * t111 + (-t113 * t39 - t114 * t38) * t141 + t38 * t50 + t40 * t44 + m(7) * (t38 * t56 - t55 * t98) - mrSges(7,3) * t98 + t22; 0.2e1 * t44 * t56 + t46 * t71 + t47 * t69 + (t118 - t120) * qJD(6); t71 * t10 + t69 * t11 + (t71 * t24 - t69 * t25) * qJD(6) + m(7) * (t1 * t69 + t2 * t71 + (-t3 * t69 + t4 * t71) * qJD(6)) - t145; m(7) * (qJD(6) * t80 + t12 * t69 + t13 * t71); 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t74; mrSges(7,1) * t13 - mrSges(7,2) * t12; t9; t58 + (t50 * t55 - t132) * qJD(6); -t44; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
