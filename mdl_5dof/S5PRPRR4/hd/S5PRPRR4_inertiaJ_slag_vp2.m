% Calculate joint inertia matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:42
% EndTime: 2019-12-05 15:49:45
% DurationCPUTime: 0.36s
% Computational Cost: add. (333->126), mult. (817->188), div. (0->0), fcn. (756->10), ass. (0->61)
t36 = sin(pkin(10));
t26 = pkin(2) * t36 + pkin(7);
t69 = 0.2e1 * t26;
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t17 = -mrSges(6,1) * t43 + mrSges(6,2) * t40;
t68 = -m(6) * pkin(4) - mrSges(5,1) + t17;
t37 = sin(pkin(5));
t38 = cos(pkin(10));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t10 = (t36 * t45 + t38 * t42) * t37;
t39 = cos(pkin(5));
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t5 = t10 * t41 - t39 * t44;
t67 = t5 ^ 2;
t8 = (t36 * t42 - t38 * t45) * t37;
t66 = t8 ^ 2;
t65 = t43 / 0.2e1;
t64 = t44 * pkin(4);
t63 = t44 * t5;
t62 = t5 * t41;
t61 = Ifges(6,4) * t40;
t60 = Ifges(6,4) * t43;
t59 = Ifges(6,6) * t44;
t58 = t26 * t44;
t57 = t40 * t41;
t56 = t41 * t43;
t55 = t40 ^ 2 + t43 ^ 2;
t33 = t41 ^ 2;
t35 = t44 ^ 2;
t54 = t33 + t35;
t27 = -pkin(2) * t38 - pkin(3);
t53 = t55 * t41;
t7 = t10 * t44 + t39 * t41;
t1 = -t40 * t7 + t43 * t8;
t2 = t40 * t8 + t43 * t7;
t52 = -t1 * t40 + t2 * t43;
t14 = -pkin(8) * t41 + t27 - t64;
t3 = t14 * t43 - t40 * t58;
t4 = t14 * t40 + t43 * t58;
t51 = -t3 * t40 + t4 * t43;
t50 = t44 * t7 + t62;
t18 = -t44 * mrSges(5,1) + t41 * mrSges(5,2);
t49 = mrSges(6,1) * t40 + mrSges(6,2) * t43;
t15 = mrSges(6,2) * t44 - mrSges(6,3) * t57;
t16 = -mrSges(6,1) * t44 - mrSges(6,3) * t56;
t48 = t43 * t15 - t40 * t16;
t31 = t39 ^ 2;
t29 = Ifges(6,5) * t40;
t28 = Ifges(6,6) * t43;
t24 = t26 ^ 2;
t23 = Ifges(6,5) * t56;
t21 = t33 * t24;
t20 = Ifges(6,1) * t40 + t60;
t19 = Ifges(6,2) * t43 + t61;
t13 = t49 * t41;
t12 = -Ifges(6,5) * t44 + (Ifges(6,1) * t43 - t61) * t41;
t11 = -t59 + (-Ifges(6,2) * t40 + t60) * t41;
t6 = [m(2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t67) + m(5) * (t7 ^ 2 + t66 + t67) + m(4) * (t10 ^ 2 + t31 + t66) + m(3) * (t31 + (t42 ^ 2 + t45 ^ 2) * t37 ^ 2); -t10 * mrSges(4,2) + t1 * t16 + t5 * t13 + t2 * t15 + (-mrSges(4,1) + t18) * t8 + (mrSges(3,1) * t45 - mrSges(3,2) * t42) * t37 + t50 * mrSges(5,3) + m(6) * (t1 * t3 + t2 * t4 + t26 * t62) + m(5) * (t26 * t50 + t27 * t8) + m(4) * (t10 * t36 - t38 * t8) * pkin(2); 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16 + 0.2e1 * t27 * t18 + Ifges(3,3) + Ifges(4,3) + (-t23 + (Ifges(6,3) + Ifges(5,2)) * t44) * t44 + m(6) * (t3 ^ 2 + t4 ^ 2 + t21) + m(5) * (t24 * t35 + t27 ^ 2 + t21) + m(4) * (t36 ^ 2 + t38 ^ 2) * pkin(2) ^ 2 + (Ifges(5,1) * t41 + 0.2e1 * Ifges(5,4) * t44 + t43 * t12 + t13 * t69 + (-t11 + t59) * t40) * t41 + 0.2e1 * (t38 * mrSges(4,1) - t36 * mrSges(4,2)) * pkin(2) + t54 * mrSges(5,3) * t69; m(4) * t39 + m(6) * (t41 * t52 - t63) + m(5) * (t41 * t7 - t63); -t44 * t13 + (m(6) * (t51 - t58) + t48) * t41; m(4) + m(5) * t54 + m(6) * (t33 * t55 + t35); -t7 * mrSges(5,2) + (m(6) * pkin(8) + mrSges(6,3)) * t52 + t68 * t5; t40 * t12 / 0.2e1 + t11 * t65 - pkin(4) * t13 + (-t26 * mrSges(5,2) - t29 / 0.2e1 - t28 / 0.2e1 + Ifges(5,6)) * t44 + t51 * mrSges(6,3) + (m(6) * t51 + t48) * pkin(8) + (t20 * t65 - t40 * t19 / 0.2e1 + Ifges(5,5) + t68 * t26) * t41; -t44 * t17 + m(6) * (pkin(8) * t53 + t64) + mrSges(6,3) * t53 - t18; Ifges(5,3) + m(6) * (pkin(8) ^ 2 * t55 + pkin(4) ^ 2) - 0.2e1 * pkin(4) * t17 + t40 * t20 + t43 * t19 + 0.2e1 * t55 * pkin(8) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2; mrSges(6,1) * t3 - mrSges(6,2) * t4 - Ifges(6,6) * t57 - Ifges(6,3) * t44 + t23; -t13; -pkin(8) * t49 + t28 + t29; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
