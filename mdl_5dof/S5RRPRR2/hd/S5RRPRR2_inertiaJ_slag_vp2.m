% Calculate joint inertia matrix for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:06
% EndTime: 2019-12-05 18:27:07
% DurationCPUTime: 0.39s
% Computational Cost: add. (1057->133), mult. (1983->199), div. (0->0), fcn. (2201->8), ass. (0->56)
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t35 = -t48 * t52 + t49 * t55;
t36 = t48 * t55 + t49 * t52;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t26 = t54 * t35 - t51 * t36;
t27 = t51 * t35 + t54 * t36;
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t12 = t53 * t26 - t50 * t27;
t76 = 0.2e1 * t12;
t75 = 0.2e1 * t26;
t74 = 0.2e1 * t35;
t73 = pkin(2) * t48;
t43 = t49 * pkin(2) + pkin(3);
t32 = t54 * t43 - t51 * t73;
t31 = pkin(4) + t32;
t33 = t51 * t43 + t54 * t73;
t20 = t50 * t31 + t53 * t33;
t72 = t20 * mrSges(6,2);
t71 = t32 * mrSges(5,1);
t70 = t33 * mrSges(5,2);
t69 = t50 * mrSges(6,2);
t68 = Ifges(5,3) + Ifges(6,3);
t67 = -qJ(3) - pkin(6);
t40 = t67 * t52;
t41 = t67 * t55;
t28 = t49 * t40 + t48 * t41;
t21 = -t36 * pkin(7) + t28;
t29 = t48 * t40 - t49 * t41;
t22 = t35 * pkin(7) + t29;
t8 = t51 * t21 + t54 * t22;
t66 = t52 ^ 2 + t55 ^ 2;
t65 = pkin(4) * t69;
t44 = -t55 * pkin(2) - pkin(1);
t13 = t50 * t26 + t53 * t27;
t64 = -t12 * mrSges(6,1) + t13 * mrSges(6,2);
t63 = -t35 * mrSges(4,1) + t36 * mrSges(4,2);
t62 = -t26 * mrSges(5,1) + t27 * mrSges(5,2);
t7 = t54 * t21 - t51 * t22;
t19 = t53 * t31 - t50 * t33;
t18 = t19 * mrSges(6,1);
t61 = Ifges(6,3) + t18 - t72;
t4 = -t27 * pkin(8) + t7;
t5 = t26 * pkin(8) + t8;
t2 = t53 * t4 - t50 * t5;
t3 = t50 * t4 + t53 * t5;
t60 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t12;
t30 = -t35 * pkin(3) + t44;
t59 = t7 * mrSges(5,1) - t8 * mrSges(5,2) + Ifges(5,5) * t27 + Ifges(5,6) * t26 + t60;
t45 = t53 * pkin(4) * mrSges(6,1);
t14 = -t26 * pkin(4) + t30;
t1 = [-0.2e1 * pkin(1) * (-t55 * mrSges(3,1) + t52 * mrSges(3,2)) + t52 * (Ifges(3,1) * t52 + Ifges(3,4) * t55) + t55 * (Ifges(3,4) * t52 + Ifges(3,2) * t55) + 0.2e1 * t44 * t63 + Ifges(4,2) * t35 ^ 2 + Ifges(5,2) * t26 ^ 2 + 0.2e1 * t30 * t62 + Ifges(6,2) * t12 ^ 2 + 0.2e1 * t14 * t64 + Ifges(2,3) + t3 * mrSges(6,3) * t76 + t8 * mrSges(5,3) * t75 + t29 * mrSges(4,3) * t74 + 0.2e1 * t66 * pkin(6) * mrSges(3,3) + (-0.2e1 * t28 * mrSges(4,3) + Ifges(4,1) * t36 + Ifges(4,4) * t74) * t36 + (-0.2e1 * t7 * mrSges(5,3) + Ifges(5,1) * t27 + Ifges(5,4) * t75) * t27 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t13 + Ifges(6,4) * t76) * t13 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t30 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2 + t44 ^ 2) + m(3) * (t66 * pkin(6) ^ 2 + pkin(1) ^ 2); t59 + (m(4) * (t28 * t49 + t29 * t48) + (t48 * t35 - t49 * t36) * mrSges(4,3)) * pkin(2) + m(6) * (t19 * t2 + t20 * t3) + m(5) * (t32 * t7 + t33 * t8) + (-t52 * mrSges(3,1) - t55 * mrSges(3,2)) * pkin(6) + (t20 * t12 - t19 * t13) * mrSges(6,3) + (t33 * t26 - t32 * t27) * mrSges(5,3) + Ifges(3,6) * t55 + Ifges(3,5) * t52 + Ifges(4,6) * t35 + Ifges(4,5) * t36 + t28 * mrSges(4,1) - t29 * mrSges(4,2); 0.2e1 * t71 - 0.2e1 * t70 - 0.2e1 * t72 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t18 + m(6) * (t19 ^ 2 + t20 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + t68 + (0.2e1 * t49 * mrSges(4,1) - 0.2e1 * t48 * mrSges(4,2) + m(4) * (t48 ^ 2 + t49 ^ 2) * pkin(2)) * pkin(2); m(4) * t44 + m(5) * t30 + m(6) * t14 + t62 + t63 + t64; 0; m(4) + m(5) + m(6); (m(6) * (t2 * t53 + t3 * t50) + (t50 * t12 - t53 * t13) * mrSges(6,3)) * pkin(4) + t59; t71 - t70 + Ifges(5,3) + t45 + (m(6) * (t19 * t53 + t20 * t50) - t69) * pkin(4) + t61; 0; -0.2e1 * t65 + 0.2e1 * t45 + m(6) * (t50 ^ 2 + t53 ^ 2) * pkin(4) ^ 2 + t68; t60; t61; 0; Ifges(6,3) + t45 - t65; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
