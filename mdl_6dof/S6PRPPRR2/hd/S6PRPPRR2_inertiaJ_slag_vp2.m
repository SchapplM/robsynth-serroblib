% Calculate joint inertia matrix for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:54:02
% EndTime: 2018-11-23 14:54:03
% DurationCPUTime: 0.46s
% Computational Cost: add. (456->153), mult. (1057->216), div. (0->0), fcn. (949->10), ass. (0->69)
t85 = m(5) + m(4);
t45 = cos(pkin(11));
t33 = -t45 * pkin(2) - pkin(3);
t29 = -pkin(8) + t33;
t84 = -0.2e1 * t29;
t82 = m(7) * pkin(9) + mrSges(7,3);
t47 = sin(qJ(6));
t50 = cos(qJ(6));
t23 = -t50 * mrSges(7,1) + t47 * mrSges(7,2);
t81 = m(7) * pkin(5) + mrSges(6,1) - t23;
t43 = sin(pkin(11));
t44 = sin(pkin(6));
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t10 = (t43 * t49 - t45 * t52) * t44;
t46 = cos(pkin(6));
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t5 = -t51 * t10 + t46 * t48;
t80 = t5 ^ 2;
t12 = (-t43 * t52 - t45 * t49) * t44;
t9 = t12 ^ 2;
t31 = t43 * pkin(2) + qJ(4);
t79 = t31 ^ 2;
t78 = t50 / 0.2e1;
t77 = t48 * pkin(5);
t76 = t48 * t5;
t75 = t51 * t5;
t40 = t48 ^ 2;
t42 = t51 ^ 2;
t64 = t42 + t40;
t74 = m(6) * t64 + m(5);
t73 = Ifges(7,4) * t47;
t72 = Ifges(7,4) * t50;
t71 = Ifges(7,6) * t48;
t70 = t29 * t48;
t69 = t31 * t12;
t68 = t51 * mrSges(7,3);
t24 = t48 * mrSges(6,1) + t51 * mrSges(6,2);
t67 = t24 + mrSges(5,3);
t66 = Ifges(7,5) * t50 * t51 + Ifges(7,3) * t48;
t65 = t47 ^ 2 + t50 ^ 2;
t61 = t65 * t51;
t60 = t64 * mrSges(6,3);
t7 = t48 * t10 + t46 * t51;
t1 = -t50 * t12 - t47 * t7;
t2 = -t47 * t12 + t50 * t7;
t59 = -t1 * t47 + t2 * t50;
t16 = -t51 * pkin(9) + t31 + t77;
t3 = t50 * t16 - t47 * t70;
t4 = t47 * t16 + t50 * t70;
t58 = -t3 * t47 + t4 * t50;
t57 = t48 * t7 - t75;
t56 = mrSges(7,1) * t47 + mrSges(7,2) * t50;
t19 = -t48 * mrSges(7,2) - t47 * t68;
t20 = t48 * mrSges(7,1) - t50 * t68;
t55 = t50 * t19 - t47 * t20;
t38 = t46 ^ 2;
t36 = Ifges(7,5) * t47;
t35 = Ifges(7,6) * t50;
t28 = t29 ^ 2;
t26 = Ifges(7,1) * t47 + t72;
t25 = Ifges(7,2) * t50 + t73;
t22 = t42 * t29;
t21 = t42 * t28;
t17 = t56 * t51;
t15 = Ifges(7,5) * t48 + (Ifges(7,1) * t50 - t73) * t51;
t14 = t71 + (-Ifges(7,2) * t47 + t72) * t51;
t6 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t80) + m(6) * (t7 ^ 2 + t80 + t9) + m(3) * (t38 + (t49 ^ 2 + t52 ^ 2) * t44 ^ 2) + t85 * (t10 ^ 2 + t38 + t9); t1 * t20 + t5 * t17 + t2 * t19 + (t52 * mrSges(3,1) - t49 * mrSges(3,2)) * t44 + (-mrSges(4,1) + mrSges(5,2)) * t10 - t57 * mrSges(6,3) - (-mrSges(4,2) + t67) * t12 + m(7) * (t3 * t1 + t4 * t2 - t29 * t75) + m(6) * (t29 * t57 - t69) + m(5) * (t33 * t10 - t69) + m(4) * (-t10 * t45 - t12 * t43) * pkin(2); 0.2e1 * t33 * mrSges(5,2) + 0.2e1 * t4 * t19 + 0.2e1 * t3 * t20 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + (Ifges(6,2) * t48 + t66) * t48 + (Ifges(6,1) * t51 - 0.2e1 * Ifges(6,4) * t48 + t50 * t15 + t17 * t84 + (-t14 - t71) * t47) * t51 + m(7) * (t3 ^ 2 + t4 ^ 2 + t21) + m(6) * (t40 * t28 + t21 + t79) + m(5) * (t33 ^ 2 + t79) + m(4) * (t43 ^ 2 + t45 ^ 2) * pkin(2) ^ 2 + 0.2e1 * t67 * t31 + 0.2e1 * (t45 * mrSges(4,1) - t43 * mrSges(4,2)) * pkin(2) + t60 * t84; t85 * t46 + m(7) * (t51 * t59 + t76) + m(6) * (t51 * t7 + t76); t48 * t17 + (m(7) * (t58 - t70) + t55) * t51; m(4) + m(7) * (t42 * t65 + t40) + t74; m(7) * (t48 * t59 - t75) + m(6) * t57 + m(5) * t10; -t51 * t17 + mrSges(5,2) + t55 * t48 - t60 + m(7) * (t48 * t58 + t22) + m(6) * (t40 * t29 + t22) + m(5) * t33; m(7) * (-0.1e1 + t65) * t51 * t48; m(7) * (t40 * t65 + t42) + t74; -t7 * mrSges(6,2) - t81 * t5 + t82 * t59; t47 * t15 / 0.2e1 + t14 * t78 - pkin(5) * t17 + (-t29 * mrSges(6,2) + t36 / 0.2e1 + t35 / 0.2e1 - Ifges(6,6)) * t48 + t58 * mrSges(7,3) + (m(7) * t58 + t55) * pkin(9) + (t26 * t78 - t47 * t25 / 0.2e1 + Ifges(6,5) + t81 * t29) * t51; t48 * t23 + m(7) * (pkin(9) * t61 - t77) + mrSges(7,3) * t61 - t24; t81 * t51 + (t65 * t82 - mrSges(6,2)) * t48; Ifges(6,3) - 0.2e1 * pkin(5) * t23 + t47 * t26 + t50 * t25 + m(7) * (pkin(9) ^ 2 * t65 + pkin(5) ^ 2) + 0.2e1 * t65 * pkin(9) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2); -Ifges(7,6) * t47 * t51 + t3 * mrSges(7,1) - t4 * mrSges(7,2) + t66; -t17; -t56 * t48; -pkin(9) * t56 + t35 + t36; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
