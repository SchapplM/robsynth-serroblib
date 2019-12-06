% Calculate joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:30
% EndTime: 2019-12-05 16:21:32
% DurationCPUTime: 0.34s
% Computational Cost: add. (451->117), mult. (974->175), div. (0->0), fcn. (983->8), ass. (0->48)
t59 = m(5) * pkin(3);
t38 = sin(pkin(9));
t39 = cos(pkin(9));
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t24 = -t38 * t41 + t39 * t44;
t25 = t38 * t44 + t39 * t41;
t15 = -t24 * mrSges(5,1) + t25 * mrSges(5,2);
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t13 = t43 * t24 - t40 * t25;
t14 = t40 * t24 + t43 * t25;
t4 = -t13 * mrSges(6,1) + t14 * mrSges(6,2);
t58 = -t15 - t4;
t36 = t44 ^ 2;
t57 = m(5) + m(6);
t56 = pkin(3) * t38;
t32 = t39 * pkin(3) + pkin(4);
t21 = t43 * t32 - t40 * t56;
t55 = t21 * mrSges(6,1);
t22 = t40 * t32 + t43 * t56;
t54 = t22 * mrSges(6,2);
t53 = -qJ(4) - pkin(6);
t29 = t53 * t41;
t31 = t53 * t44;
t17 = t38 * t29 - t39 * t31;
t52 = t41 ^ 2 + t36;
t42 = sin(qJ(2));
t19 = t25 * t42;
t20 = t24 * t42;
t6 = -t43 * t19 - t40 * t20;
t7 = -t40 * t19 + t43 * t20;
t51 = t6 * mrSges(6,1) - t7 * mrSges(6,2);
t33 = -t44 * pkin(3) - pkin(2);
t50 = t52 * mrSges(4,3);
t16 = t39 * t29 + t38 * t31;
t49 = -t41 * mrSges(4,1) - t44 * mrSges(4,2);
t8 = -t25 * pkin(7) + t16;
t9 = t24 * pkin(7) + t17;
t2 = -t40 * t9 + t43 * t8;
t3 = t40 * t8 + t43 * t9;
t48 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t14 + Ifges(6,6) * t13;
t45 = cos(qJ(2));
t37 = t45 ^ 2;
t35 = t42 ^ 2;
t30 = -t44 * mrSges(4,1) + t41 * mrSges(4,2);
t18 = -t24 * pkin(4) + t33;
t1 = [m(2) + m(6) * (t6 ^ 2 + t7 ^ 2 + t37) + m(5) * (t19 ^ 2 + t20 ^ 2 + t37) + m(4) * (t52 * t35 + t37) + m(3) * (t35 + t37); (t7 * t13 - t6 * t14) * mrSges(6,3) + (t19 * t25 + t20 * t24) * mrSges(5,3) + (-mrSges(3,2) + t50) * t42 + (mrSges(3,1) - t30 + t58) * t45 + m(6) * (-t18 * t45 + t2 * t6 + t3 * t7) + m(5) * (-t16 * t19 + t17 * t20 - t33 * t45) + m(4) * (t52 * t42 * pkin(6) + t45 * pkin(2)); Ifges(4,2) * t36 - 0.2e1 * pkin(2) * t30 + 0.2e1 * t33 * t15 + 0.2e1 * t18 * t4 + Ifges(3,3) + (Ifges(4,1) * t41 + 0.2e1 * Ifges(4,4) * t44) * t41 + (-0.2e1 * t16 * mrSges(5,3) + Ifges(5,1) * t25) * t25 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t14) * t14 + 0.2e1 * pkin(6) * t50 + (0.2e1 * t17 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t25 + Ifges(5,2) * t24) * t24 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t14 + Ifges(6,2) * t13) * t13 + m(6) * (t18 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t33 ^ 2) + m(4) * (t52 * pkin(6) ^ 2 + pkin(2) ^ 2); -t19 * mrSges(5,1) - t20 * mrSges(5,2) + t49 * t42 + m(6) * (t21 * t6 + t22 * t7) + (-t19 * t39 + t20 * t38) * t59 + t51; m(6) * (t21 * t2 + t22 * t3) + t16 * mrSges(5,1) - t17 * mrSges(5,2) + Ifges(5,5) * t25 + Ifges(5,6) * t24 + Ifges(4,5) * t41 + Ifges(4,6) * t44 + t49 * pkin(6) + (t22 * t13 - t21 * t14) * mrSges(6,3) + (m(5) * (t16 * t39 + t17 * t38) + (t38 * t24 - t39 * t25) * mrSges(5,3)) * pkin(3) + t48; 0.2e1 * t55 - 0.2e1 * t54 + Ifges(4,3) + Ifges(5,3) + Ifges(6,3) + m(6) * (t21 ^ 2 + t22 ^ 2) + (0.2e1 * t39 * mrSges(5,1) - 0.2e1 * t38 * mrSges(5,2) + (t38 ^ 2 + t39 ^ 2) * t59) * pkin(3); -t57 * t45; m(5) * t33 + m(6) * t18 - t58; 0; t57; t51; t48; Ifges(6,3) - t54 + t55; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
