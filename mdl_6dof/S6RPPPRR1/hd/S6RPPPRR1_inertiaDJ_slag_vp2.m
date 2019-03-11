% Calculate time derivative of joint inertia matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:51
% EndTime: 2019-03-09 01:29:53
% DurationCPUTime: 0.84s
% Computational Cost: add. (656->174), mult. (1348->269), div. (0->0), fcn. (884->6), ass. (0->83)
t39 = sin(qJ(6));
t41 = cos(qJ(6));
t70 = t39 ^ 2 + t41 ^ 2;
t92 = m(7) * (-0.1e1 + t70);
t22 = -mrSges(7,1) * t41 + mrSges(7,2) * t39;
t90 = -mrSges(6,1) + t22;
t40 = sin(qJ(5));
t34 = t40 ^ 2;
t42 = cos(qJ(5));
t36 = t42 ^ 2;
t89 = t34 - t36;
t88 = 2 * qJD(4);
t64 = qJD(5) * t42;
t87 = t40 * t64 * t92;
t86 = -t39 / 0.2e1;
t85 = pkin(5) * t42;
t84 = pkin(8) * t40;
t83 = pkin(8) * t42;
t82 = t40 * pkin(5);
t81 = Ifges(7,4) * t39;
t80 = Ifges(7,4) * t41;
t79 = Ifges(7,5) * t39;
t78 = Ifges(7,6) * t39;
t77 = Ifges(7,6) * t40;
t76 = Ifges(7,6) * t41;
t75 = t39 * t40;
t74 = t40 * t41;
t72 = t42 * mrSges(7,3);
t19 = -t40 * mrSges(7,2) - t39 * t72;
t73 = t41 * t19;
t65 = qJD(5) * t40;
t56 = t39 * t65;
t71 = Ifges(7,6) * t56 + Ifges(7,3) * t64;
t68 = qJD(3) * t34;
t67 = qJD(3) * t40;
t27 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t66 = qJD(4) * t27;
t63 = qJD(6) * t39;
t62 = qJD(6) * t40;
t61 = qJD(6) * t41;
t60 = qJD(6) * t42;
t32 = t36 * qJD(3);
t28 = sin(pkin(9)) * pkin(1) + qJ(3);
t26 = -pkin(7) + t28;
t59 = t26 * t64;
t58 = t39 * t64;
t57 = t41 * t60;
t55 = t41 * t64;
t54 = -Ifges(7,5) * t41 + (2 * Ifges(6,4));
t13 = t27 + t82 - t83;
t45 = qJD(6) * t13 + t59 + t67;
t51 = -t26 * t62 + qJD(4) + (t84 + t85) * qJD(5);
t1 = t39 * t51 + t41 * t45;
t6 = t13 * t41 - t26 * t75;
t53 = -qJD(6) * t6 + t1;
t52 = m(7) * pkin(5) - t90;
t50 = mrSges(7,1) * t39 + mrSges(7,2) * t41;
t49 = Ifges(7,1) * t41 - t81;
t24 = Ifges(7,1) * t39 + t80;
t48 = -Ifges(7,2) * t39 + t80;
t23 = Ifges(7,2) * t41 + t81;
t47 = t39 * t60 + t41 * t65;
t46 = t56 - t57;
t2 = -t39 * t45 + t41 * t51;
t7 = t13 * t39 + t26 * t74;
t44 = t1 * t41 - t2 * t39 - t6 * t61 - t63 * t7;
t10 = -mrSges(7,2) * t64 + mrSges(7,3) * t46;
t14 = t50 * t42;
t20 = t40 * mrSges(7,1) - t41 * t72;
t9 = mrSges(7,1) * t64 + mrSges(7,3) * t47;
t43 = qJD(5) * t14 + t41 * t10 - t19 * t63 - t20 * t61 - t39 * t9;
t31 = Ifges(7,5) * t61;
t29 = mrSges(6,2) * t65;
t21 = t26 * t32;
t18 = t49 * qJD(6);
t17 = t48 * qJD(6);
t16 = t50 * qJD(6);
t12 = Ifges(7,5) * t40 + t42 * t49;
t11 = t42 * t48 + t77;
t5 = -mrSges(7,1) * t46 - mrSges(7,2) * t47;
t4 = -t24 * t60 + (Ifges(7,5) * t42 - t40 * t49) * qJD(5);
t3 = -t23 * t60 + (Ifges(7,6) * t42 - t48 * t40) * qJD(5);
t8 = [(mrSges(5,3) * t88) + 0.2e1 * t1 * t19 + 0.2e1 * t7 * t10 + 0.2e1 * t2 * t20 - 0.2e1 * t27 * t29 + 0.2e1 * t6 * t9 + 0.2e1 * (m(4) * t28 + mrSges(5,2) + mrSges(4,3) + (-t34 - t36) * mrSges(6,3)) * qJD(3) + 0.2e1 * m(7) * (t7 * t1 + t6 * t2 + t21) + 0.2e1 * m(5) * (qJD(3) * t28 + t66) + 0.2e1 * m(6) * (t26 * t68 + t21 + t66) + (mrSges(6,1) * t88 + (t39 * t11 - t41 * t12 + 0.2e1 * t26 * t14 + t40 * t54) * qJD(5) + t71) * t40 + (mrSges(6,2) * t88 - 0.2e1 * qJD(3) * t14 - 0.2e1 * t26 * t5 - t39 * t3 + t41 * t4 + (-t39 * t12 + t40 * (-t76 - t79) - t41 * t11) * qJD(6) + (0.2e1 * t27 * mrSges(6,1) + (-t54 - t78) * t42 + (-0.2e1 * m(7) * t26 ^ 2 - (2 * Ifges(6,1)) + (2 * Ifges(6,2)) + Ifges(7,3)) * t40) * qJD(5)) * t42; t40 * t5 + (t20 * t75 - t40 * t73 + m(7) * (t89 * t26 + t6 * t75 - t7 * t74)) * qJD(5) + (m(7) * (t44 - t67) + t43) * t42; -0.2e1 * t87; m(7) * (-t1 * t39 - t2 * t41 + (t39 * t6 - t41 * t7) * qJD(6)) - t19 * t61 - t39 * t10 + t20 * t63 - t41 * t9 - mrSges(6,1) * t64 + t29 + (-m(6) - m(5)) * qJD(4); 0; 0; m(5) * qJD(3) + (-t5 + (-t39 * t20 + t73) * qJD(5)) * t42 + m(7) * (t7 * t55 - t6 * t58 + t32) + m(6) * (t32 + t68) + (m(7) * (t44 - 0.2e1 * t59) + t43) * t40; -t89 * qJD(5) * t92; 0; 0.2e1 * t87; -pkin(5) * t5 + (t31 / 0.2e1 - qJD(3) * mrSges(6,2) + (-t26 * t52 - Ifges(6,5)) * qJD(5)) * t40 + (t4 / 0.2e1 - t2 * mrSges(7,3) + t23 * t65 / 0.2e1 + (-t77 / 0.2e1 - t11 / 0.2e1 - t7 * mrSges(7,3)) * qJD(6) + (m(7) * (-qJD(6) * t7 - t2) - qJD(6) * t19 - t9) * pkin(8)) * t39 + (t3 / 0.2e1 + qJD(6) * t12 / 0.2e1 - t24 * t65 / 0.2e1 + t53 * mrSges(7,3) + (m(7) * t53 - qJD(6) * t20 + t10) * pkin(8)) * t41 + (t41 * t18 / 0.2e1 + t17 * t86 - t26 * t16 + (-t41 * t23 / 0.2e1 + t24 * t86) * qJD(6) + (t79 / 0.2e1 + t76 / 0.2e1 - Ifges(6,6) - t26 * mrSges(6,2)) * qJD(5) + t52 * qJD(3)) * t42; t40 * t16 + t29 + (m(7) * (-t70 * t84 - t85) - t70 * t40 * mrSges(7,3) + t90 * t42) * qJD(5); 0; -t42 * t16 + (-t42 * mrSges(6,2) + m(7) * (t70 * t83 - t82) + t70 * t72 + t90 * t40) * qJD(5); -0.2e1 * pkin(5) * t16 + t17 * t41 + t18 * t39 + (-t23 * t39 + t24 * t41) * qJD(6); t2 * mrSges(7,1) - t1 * mrSges(7,2) - Ifges(7,5) * t47 - Ifges(7,6) * t57 + t71; -t5; t16; (t39 * t62 - t55) * mrSges(7,2) + (-t40 * t61 - t58) * mrSges(7,1); t31 + (pkin(8) * t22 - t78) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
