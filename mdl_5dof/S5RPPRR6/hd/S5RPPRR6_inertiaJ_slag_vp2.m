% Calculate joint inertia matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:42
% DurationCPUTime: 0.32s
% Computational Cost: add. (461->111), mult. (873->164), div. (0->0), fcn. (833->8), ass. (0->55)
t38 = cos(pkin(9));
t36 = sin(pkin(9));
t41 = sin(qJ(4));
t56 = t41 * t36;
t62 = cos(qJ(4));
t18 = -t62 * t38 + t56;
t50 = t62 * t36;
t20 = t41 * t38 + t50;
t40 = sin(qJ(5));
t59 = t20 * t40;
t10 = -t18 * mrSges(6,2) - mrSges(6,3) * t59;
t42 = cos(qJ(5));
t58 = t20 * t42;
t11 = t18 * mrSges(6,1) - mrSges(6,3) * t58;
t71 = t42 * t10 - t40 * t11;
t37 = sin(pkin(8));
t26 = t37 * pkin(1) + qJ(3);
t63 = pkin(6) + t26;
t14 = t63 * t38;
t6 = t41 * t14 + t63 * t50;
t70 = t6 ^ 2;
t69 = 0.2e1 * t6;
t68 = t18 ^ 2;
t33 = t38 ^ 2;
t39 = cos(pkin(8));
t28 = -t39 * pkin(1) - pkin(2);
t21 = -t38 * pkin(3) + t28;
t67 = 0.2e1 * t21;
t66 = t42 / 0.2e1;
t65 = pkin(4) * t18;
t64 = t6 * t18;
t61 = Ifges(6,4) * t40;
t60 = Ifges(6,4) * t42;
t54 = Ifges(6,5) * t58 + Ifges(6,3) * t18;
t53 = Ifges(6,5) * t40 + Ifges(6,6) * t42;
t52 = t36 ^ 2 + t33;
t51 = t40 ^ 2 + t42 ^ 2;
t49 = t51 * t20;
t48 = -t38 * mrSges(4,1) + t36 * mrSges(4,2);
t15 = t20 * mrSges(5,2);
t47 = -t18 * mrSges(5,1) - t15;
t5 = -t20 * pkin(7) + t21 + t65;
t8 = t62 * t14 - t63 * t56;
t1 = -t40 * t8 + t42 * t5;
t2 = t40 * t5 + t42 * t8;
t46 = -t1 * t40 + t2 * t42;
t45 = mrSges(6,1) * t40 + mrSges(6,2) * t42;
t24 = Ifges(6,1) * t40 + t60;
t23 = Ifges(6,2) * t42 + t61;
t22 = -t42 * mrSges(6,1) + t40 * mrSges(6,2);
t17 = t20 ^ 2;
t9 = t45 * t20;
t4 = Ifges(6,5) * t18 + (Ifges(6,1) * t42 - t61) * t20;
t3 = Ifges(6,6) * t18 + (-Ifges(6,2) * t40 + t60) * t20;
t7 = [Ifges(2,3) + Ifges(3,3) + t9 * t69 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + t15 * t67 + 0.2e1 * t28 * t48 + Ifges(4,2) * t33 + (Ifges(4,1) * t36 + 0.2e1 * Ifges(4,4) * t38) * t36 + (mrSges(5,1) * t67 - 0.2e1 * t8 * mrSges(5,3) + Ifges(5,2) * t18 + t54) * t18 + (mrSges(5,3) * t69 + Ifges(5,1) * t20 - t40 * t3 + t42 * t4 + (-Ifges(6,6) * t40 - (2 * Ifges(5,4))) * t18) * t20 + m(6) * (t1 ^ 2 + t2 ^ 2 + t70) + m(5) * (t21 ^ 2 + t8 ^ 2 + t70) + m(4) * (t52 * t26 ^ 2 + t28 ^ 2) + m(3) * (t37 ^ 2 + t39 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t39 * mrSges(3,1) - t37 * mrSges(3,2)) * pkin(1) + 0.2e1 * t52 * t26 * mrSges(4,3); t18 * t9 + t71 * t20 + m(5) * (t8 * t20 + t64) + m(6) * (t46 * t20 + t64); m(3) + m(5) * (t17 + t68) + m(6) * (t51 * t17 + t68) + m(4) * t52; t40 * t10 + t42 * t11 + m(6) * (t42 * t1 + t40 * t2) + m(5) * t21 + m(4) * t28 - t47 + t48; 0; m(6) * t51 + m(4) + m(5); t40 * t4 / 0.2e1 + t3 * t66 - pkin(4) * t9 - t8 * mrSges(5,2) + t46 * mrSges(6,3) + (-t40 * t23 / 0.2e1 + t24 * t66 + Ifges(5,5)) * t20 + (t53 / 0.2e1 - Ifges(5,6)) * t18 + (-m(6) * pkin(4) - mrSges(5,1) + t22) * t6 + (m(6) * t46 + t71) * pkin(7); m(6) * (pkin(7) * t49 - t65) + t18 * t22 + mrSges(6,3) * t49 + t47; 0; Ifges(5,3) + t40 * t24 + t42 * t23 + m(6) * (t51 * pkin(7) ^ 2 + pkin(4) ^ 2) - 0.2e1 * pkin(4) * t22 + 0.2e1 * t51 * pkin(7) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,6) * t59 + t54; -t9; -t22; -t45 * pkin(7) + t53; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
