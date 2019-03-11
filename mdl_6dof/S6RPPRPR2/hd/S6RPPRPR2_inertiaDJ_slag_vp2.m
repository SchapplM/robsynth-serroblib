% Calculate time derivative of joint inertia matrix for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:50
% EndTime: 2019-03-09 01:41:51
% DurationCPUTime: 0.99s
% Computational Cost: add. (1652->196), mult. (3482->286), div. (0->0), fcn. (3176->8), ass. (0->94)
t82 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t122 = m(6) + m(5);
t107 = cos(qJ(4));
t57 = sin(pkin(10));
t58 = cos(pkin(10));
t61 = sin(qJ(4));
t41 = t107 * t57 + t61 * t58;
t37 = t41 * qJD(4);
t60 = sin(qJ(6));
t100 = t60 * t37;
t101 = t57 * t61;
t87 = t107 * t58;
t40 = -t87 + t101;
t62 = cos(qJ(6));
t91 = qJD(6) * t62;
t70 = t40 * t91 + t100;
t90 = 0.2e1 * t37;
t85 = qJD(4) * t107;
t36 = qJD(4) * t101 - t58 * t85;
t73 = qJ(5) * t36 - t41 * qJD(5);
t18 = pkin(4) * t37 + t73;
t33 = t36 * mrSges(6,3);
t34 = t36 * mrSges(5,2);
t119 = -m(6) * t18 - t33 + t34;
t112 = pkin(4) + pkin(8);
t10 = t112 * t37 + t73;
t72 = -cos(pkin(9)) * pkin(1) - pkin(3) * t58 - pkin(2);
t67 = -qJ(5) * t41 + t72;
t15 = t112 * t40 + t67;
t50 = sin(pkin(9)) * pkin(1) + qJ(3);
t108 = pkin(7) + t50;
t38 = t108 * t57;
t39 = t108 * t58;
t22 = t107 * t38 + t39 * t61;
t16 = pkin(5) * t41 + t22;
t3 = -t15 * t60 + t16 * t62;
t23 = t107 * t39 - t61 * t38;
t14 = t41 * qJD(3) + t23 * qJD(4);
t9 = -t36 * pkin(5) + t14;
t1 = t3 * qJD(6) + t10 * t62 + t60 * t9;
t4 = t15 * t62 + t16 * t60;
t2 = -t4 * qJD(6) - t10 * t60 + t62 * t9;
t118 = t1 * t60 + t2 * t62;
t13 = (qJD(3) * t57 + qJD(4) * t39) * t61 - qJD(3) * t87 + t38 * t85;
t117 = -0.2e1 * t18;
t116 = -t60 / 0.2e1;
t115 = t60 / 0.2e1;
t114 = -t62 / 0.2e1;
t113 = t62 / 0.2e1;
t106 = mrSges(7,3) * t40;
t105 = Ifges(7,4) * t60;
t104 = Ifges(7,4) * t62;
t102 = t41 * Ifges(7,6);
t27 = t41 * t36;
t47 = Ifges(7,1) * t62 - t105;
t99 = t60 * t47;
t98 = t62 * t37;
t46 = -Ifges(7,2) * t60 + t104;
t97 = t62 * t46;
t96 = mrSges(6,2) - mrSges(5,1);
t93 = t60 ^ 2 + t62 ^ 2;
t92 = qJD(6) * t60;
t89 = t40 * t92;
t86 = t93 * t37;
t83 = Ifges(7,5) * t70 + Ifges(7,6) * t98 - Ifges(7,3) * t36;
t81 = t3 * t62 + t4 * t60;
t80 = t3 * t60 - t4 * t62;
t45 = mrSges(7,1) * t60 + mrSges(7,2) * t62;
t79 = Ifges(7,1) * t60 + t104;
t78 = Ifges(7,2) * t62 + t105;
t77 = -Ifges(7,5) * t60 - Ifges(7,6) * t62;
t76 = -t13 * t23 + t14 * t22;
t25 = mrSges(7,1) * t41 - t60 * t106;
t26 = -mrSges(7,2) * t41 + t62 * t106;
t75 = t62 * t25 + t60 * t26;
t74 = -t60 * t25 + t62 * t26;
t69 = t89 - t98;
t65 = -t80 * qJD(6) + t118;
t11 = mrSges(7,2) * t36 - t69 * mrSges(7,3);
t12 = -mrSges(7,1) * t36 - t70 * mrSges(7,3);
t64 = t74 * qJD(6) + t60 * t11 + t62 * t12;
t44 = t79 * qJD(6);
t43 = t78 * qJD(6);
t42 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t24 = (-mrSges(7,1) * t62 + mrSges(7,2) * t60) * t40;
t21 = pkin(4) * t40 + t67;
t20 = t41 * Ifges(7,5) + t79 * t40;
t19 = t78 * t40 + t102;
t17 = -t40 * pkin(5) + t23;
t8 = -pkin(5) * t37 - t13;
t7 = t69 * mrSges(7,1) + t70 * mrSges(7,2);
t6 = t70 * Ifges(7,1) - t69 * Ifges(7,4) - Ifges(7,5) * t36;
t5 = t70 * Ifges(7,4) - t69 * Ifges(7,2) - Ifges(7,6) * t36;
t28 = [t19 * t98 + t20 * t100 + 0.2e1 * t72 * (t37 * mrSges(5,1) - t34) + 0.2e1 * t21 * (-t37 * mrSges(6,2) + t33) + 0.2e1 * t8 * t24 + 0.2e1 * t2 * t25 + 0.2e1 * t1 * t26 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t12 + 0.2e1 * t17 * t7 + 0.2e1 * m(7) * (t1 * t4 + t17 * t8 + t2 * t3) + 0.2e1 * m(5) * t76 + 0.2e1 * m(6) * (t18 * t21 + t76) + (mrSges(6,3) * t117 + (-Ifges(5,4) - Ifges(6,6)) * t90 + t14 * t82 + (-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t36 + t83) * t41 + (mrSges(6,2) * t117 + t62 * t5 + t60 * t6 + (Ifges(5,2) + Ifges(6,3)) * t90 + t13 * t82 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) + t77) * t36 + (t62 * t20 + (-t19 - t102) * t60) * qJD(6)) * t40 + 0.2e1 * (m(4) * t50 + mrSges(4,3)) * qJD(3) * (t57 ^ 2 + t58 ^ 2) + (-t22 * t36 - t23 * t37) * t82; -t36 * t24 + t41 * t7 + t75 * t37 + t64 * t40 + m(7) * (-t17 * t36 + t81 * t37 + t65 * t40 + t41 * t8) + t122 * (-t13 * t41 + t14 * t40 + t22 * t37 - t23 * t36); 0.2e1 * m(7) * (t40 * t86 - t27) + 0.2e1 * t122 * (t40 * t37 - t27); t62 * t11 - t60 * t12 - t96 * t37 - t75 * qJD(6) + m(7) * (-t81 * qJD(6) + t1 * t62 - t2 * t60) - t119; 0; 0; qJ(5) * t7 + qJD(5) * t24 - t17 * t42 + t8 * t45 + t96 * t14 + (-mrSges(6,3) + mrSges(5,2)) * t13 + (-t2 * mrSges(7,3) - t112 * t12 + t6 / 0.2e1) * t62 + (-t112 * t11 - t1 * mrSges(7,3) - t5 / 0.2e1) * t60 + m(7) * (qJ(5) * t8 + qJD(5) * t17 - t112 * t118) + m(6) * (-pkin(4) * t14 - qJ(5) * t13 + qJD(5) * t23) + (-qJD(5) * mrSges(6,1) - t43 * t113 - t44 * t115) * t40 + (pkin(4) * mrSges(6,1) + Ifges(7,5) * t114 + Ifges(7,6) * t115 + Ifges(6,4) - Ifges(5,5)) * t36 + (t97 / 0.2e1 + t99 / 0.2e1 - qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6)) * t37 + (t19 * t114 + t20 * t116 + t41 * t77 / 0.2e1 + (t47 * t113 + t46 * t116) * t40 + t80 * mrSges(7,3) - (-m(7) * t80 + t74) * t112) * qJD(6); -t36 * t45 - t41 * t42 + (-t93 * mrSges(7,3) + t96) * t37 + m(7) * (-t112 * t86 - t73) + t119; 0; -0.2e1 * qJ(5) * t42 + t43 * t60 - t44 * t62 + (-t97 - t99) * qJD(6) + 0.2e1 * (mrSges(6,3) + t45 + (m(6) + m(7)) * qJ(5)) * qJD(5); m(6) * t14 + m(7) * t65 - t36 * mrSges(6,1) + t64; (m(6) / 0.2e1 + m(7) * t93 / 0.2e1) * t90; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t89 + t83; -t7; t42; ((mrSges(7,2) * t112 - Ifges(7,6)) * t62 + (mrSges(7,1) * t112 - Ifges(7,5)) * t60) * qJD(6); -t45 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t28(1) t28(2) t28(4) t28(7) t28(11) t28(16); t28(2) t28(3) t28(5) t28(8) t28(12) t28(17); t28(4) t28(5) t28(6) t28(9) t28(13) t28(18); t28(7) t28(8) t28(9) t28(10) t28(14) t28(19); t28(11) t28(12) t28(13) t28(14) t28(15) t28(20); t28(16) t28(17) t28(18) t28(19) t28(20) t28(21);];
Mq  = res;
