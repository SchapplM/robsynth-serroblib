% Calculate joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:27
% EndTime: 2019-12-05 16:02:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (271->121), mult. (622->176), div. (0->0), fcn. (497->8), ass. (0->59)
t67 = m(6) * pkin(8) + mrSges(6,3);
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t15 = -mrSges(6,1) * t38 + mrSges(6,2) * t35;
t66 = m(6) * pkin(4) + mrSges(5,1) - t15;
t34 = cos(pkin(5));
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t33 = sin(pkin(5));
t40 = cos(qJ(2));
t58 = t33 * t40;
t8 = t34 * t36 + t39 * t58;
t65 = t8 ^ 2;
t64 = t38 / 0.2e1;
t63 = t39 * t8;
t62 = Ifges(6,4) * t35;
t61 = Ifges(6,4) * t38;
t60 = Ifges(6,6) * t36;
t37 = sin(qJ(2));
t59 = t33 * t37;
t57 = t35 * t39;
t41 = -pkin(2) - pkin(7);
t56 = t36 * t41;
t55 = t38 * t39;
t54 = mrSges(5,1) * t36 + mrSges(5,2) * t39 + mrSges(4,3);
t53 = Ifges(6,5) * t55 + Ifges(6,3) * t36;
t52 = t35 ^ 2 + t38 ^ 2;
t29 = t36 ^ 2;
t31 = t39 ^ 2;
t51 = -t31 - t29;
t49 = t51 * mrSges(5,3);
t10 = t34 * t39 - t36 * t58;
t1 = -t10 * t35 + t38 * t59;
t2 = t10 * t38 + t35 * t59;
t48 = -t1 * t35 + t2 * t38;
t14 = pkin(4) * t36 - pkin(8) * t39 + qJ(3);
t3 = t14 * t38 - t35 * t56;
t4 = t14 * t35 + t38 * t56;
t47 = -t3 * t35 + t38 * t4;
t46 = t10 * t36 - t63;
t45 = mrSges(6,1) * t35 + mrSges(6,2) * t38;
t12 = -mrSges(6,2) * t36 - mrSges(6,3) * t57;
t13 = mrSges(6,1) * t36 - mrSges(6,3) * t55;
t44 = t38 * t12 - t35 * t13;
t42 = qJ(3) ^ 2;
t32 = t41 ^ 2;
t27 = t33 ^ 2;
t26 = Ifges(6,5) * t35;
t25 = Ifges(6,6) * t38;
t23 = t31 * t41;
t22 = t31 * t32;
t21 = t27 * t37 ^ 2;
t19 = qJ(3) * t59;
t18 = Ifges(6,1) * t35 + t61;
t17 = Ifges(6,2) * t38 + t62;
t11 = t45 * t39;
t6 = Ifges(6,5) * t36 + (Ifges(6,1) * t38 - t62) * t39;
t5 = t60 + (-Ifges(6,2) * t35 + t61) * t39;
t7 = [m(2) + m(5) * (t10 ^ 2 + t21 + t65) + m(6) * (t1 ^ 2 + t2 ^ 2 + t65) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t27 * t40 ^ 2 + t34 ^ 2 + t21); t1 * t13 + t8 * t11 + t2 * t12 - t46 * mrSges(5,3) + ((mrSges(3,1) - mrSges(4,2)) * t40 + (-mrSges(3,2) + t54) * t37) * t33 + m(4) * (pkin(2) * t58 + t19) + m(5) * (t41 * t46 + t19) + m(6) * (t1 * t3 + t2 * t4 - t41 * t63); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t4 * t12 + 0.2e1 * t3 * t13 + Ifges(4,1) + Ifges(3,3) + (Ifges(5,2) * t36 + t53) * t36 + m(6) * (t3 ^ 2 + t4 ^ 2 + t22) + m(5) * (t29 * t32 + t22 + t42) + m(4) * (pkin(2) ^ 2 + t42) + (Ifges(5,1) * t39 - 0.2e1 * Ifges(5,4) * t36 - 0.2e1 * t41 * t11 + t38 * t6 + (-t5 - t60) * t35) * t39 + 0.2e1 * t54 * qJ(3) + 0.2e1 * t41 * t49; -m(4) * t58 + m(5) * t46 + m(6) * (t36 * t48 - t63); -m(4) * pkin(2) - t39 * t11 + mrSges(4,2) + t44 * t36 + t49 + m(6) * (t36 * t47 + t23) + m(5) * (t29 * t41 + t23); m(4) - m(5) * t51 + m(6) * (t29 * t52 + t31); -t10 * mrSges(5,2) + t48 * t67 - t66 * t8; t35 * t6 / 0.2e1 + t5 * t64 - pkin(4) * t11 + (-t41 * mrSges(5,2) + t26 / 0.2e1 + t25 / 0.2e1 - Ifges(5,6)) * t36 + t47 * mrSges(6,3) + (m(6) * t47 + t44) * pkin(8) + (t18 * t64 - t35 * t17 / 0.2e1 + Ifges(5,5) + t66 * t41) * t39; t66 * t39 + (t52 * t67 - mrSges(5,2)) * t36; Ifges(5,3) - 0.2e1 * pkin(4) * t15 + m(6) * (pkin(8) ^ 2 * t52 + pkin(4) ^ 2) + t35 * t18 + t38 * t17 + 0.2e1 * t52 * pkin(8) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2; mrSges(6,1) * t3 - mrSges(6,2) * t4 - Ifges(6,6) * t57 + t53; -t45 * t36; -pkin(8) * t45 + t25 + t26; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
