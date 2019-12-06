% Calculate joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:46
% EndTime: 2019-12-05 16:47:48
% DurationCPUTime: 0.31s
% Computational Cost: add. (321->109), mult. (698->146), div. (0->0), fcn. (620->6), ass. (0->43)
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t53 = -mrSges(5,2) - mrSges(6,2);
t60 = (t39 * mrSges(5,1) + t53 * t36) * pkin(3);
t40 = cos(qJ(3));
t34 = t40 ^ 2;
t59 = 2 * mrSges(6,1);
t58 = m(6) * pkin(4);
t57 = -pkin(7) - pkin(6);
t56 = t39 * pkin(3);
t37 = sin(qJ(3));
t21 = -t36 * t37 + t39 * t40;
t55 = t36 * t21;
t52 = Ifges(5,3) + Ifges(6,3);
t27 = t57 * t37;
t28 = t57 * t40;
t9 = t36 * t27 - t39 * t28;
t51 = t37 ^ 2 + t34;
t30 = -pkin(3) * t40 - pkin(2);
t49 = t51 * mrSges(4,3);
t22 = t36 * t40 + t37 * t39;
t5 = -t21 * mrSges(6,1) + t22 * mrSges(6,2);
t8 = t39 * t27 + t28 * t36;
t3 = -qJ(5) * t22 + t8;
t48 = m(6) * t3 - t22 * mrSges(6,3);
t47 = -mrSges(4,1) * t37 - mrSges(4,2) * t40;
t38 = sin(qJ(2));
t14 = t22 * t38;
t15 = t21 * t38;
t45 = t53 * t15 + (-mrSges(5,1) - mrSges(6,1)) * t14;
t4 = qJ(5) * t21 + t9;
t44 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t9 * mrSges(5,2) - t4 * mrSges(6,2) + (Ifges(5,5) + Ifges(6,5)) * t22 + (Ifges(5,6) + Ifges(6,6)) * t21;
t43 = pkin(3) ^ 2;
t41 = cos(qJ(2));
t35 = t41 ^ 2;
t33 = t38 ^ 2;
t31 = t36 ^ 2 * t43;
t29 = pkin(4) + t56;
t26 = -mrSges(4,1) * t40 + mrSges(4,2) * t37;
t11 = -pkin(4) * t21 + t30;
t10 = t36 * pkin(3) * t15;
t6 = -mrSges(5,1) * t21 + mrSges(5,2) * t22;
t1 = [m(2) + m(3) * (t33 + t35) + m(4) * (t33 * t51 + t35) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t14 ^ 2 + t15 ^ 2 + t35); (-mrSges(3,2) + t49) * t38 + (mrSges(3,1) - t26 - t5 - t6) * t41 + m(4) * (pkin(6) * t38 * t51 + t41 * pkin(2)) + m(5) * (-t14 * t8 + t15 * t9 - t30 * t41) + m(6) * (-t11 * t41 - t14 * t3 + t15 * t4) + (mrSges(6,3) + mrSges(5,3)) * (t14 * t22 + t15 * t21); Ifges(4,2) * t34 - 0.2e1 * pkin(2) * t26 + 0.2e1 * t11 * t5 + 0.2e1 * t30 * t6 + Ifges(3,3) + (Ifges(4,1) * t37 + 0.2e1 * Ifges(4,4) * t40) * t37 + m(6) * (t11 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t30 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(4) * (pkin(6) ^ 2 * t51 + pkin(2) ^ 2) + (-0.2e1 * t8 * mrSges(5,3) - 0.2e1 * t3 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1)) * t22) * t22 + 0.2e1 * pkin(6) * t49 + (0.2e1 * t9 * mrSges(5,3) + 0.2e1 * t4 * mrSges(6,3) + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t22 + (Ifges(5,2) + Ifges(6,2)) * t21) * t21; t47 * t38 + m(5) * (-t14 * t56 + t10) + m(6) * (-t14 * t29 + t10) + t45; Ifges(4,5) * t37 + Ifges(4,6) * t40 + t48 * t29 + t47 * pkin(6) + (mrSges(6,3) * t55 + (-t22 * t39 + t55) * mrSges(5,3) + m(6) * t36 * t4 + m(5) * (t36 * t9 + t39 * t8)) * pkin(3) + t44; t29 * t59 + Ifges(4,3) + m(6) * (t29 ^ 2 + t31) + m(5) * (t39 ^ 2 * t43 + t31) + 0.2e1 * t60 + t52; -t14 * t58 + t45; pkin(4) * t48 + t44; t29 * t58 + (pkin(4) + t29) * mrSges(6,1) + t60 + t52; (t59 + t58) * pkin(4) + t52; -m(6) * t41; m(6) * t11 + t5; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
