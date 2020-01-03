% Calculate joint inertia matrix for
% S5RRPRR3
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:22
% EndTime: 2020-01-03 12:00:23
% DurationCPUTime: 0.26s
% Computational Cost: add. (432->88), mult. (770->122), div. (0->0), fcn. (566->8), ass. (0->55)
t43 = cos(qJ(5));
t37 = t43 ^ 2;
t40 = sin(qJ(5));
t26 = -t43 * mrSges(6,1) + t40 * mrSges(6,2);
t45 = cos(qJ(2));
t31 = t45 * pkin(1) + pkin(2);
t38 = sin(pkin(9));
t39 = cos(pkin(9));
t42 = sin(qJ(2));
t66 = pkin(1) * t42;
t20 = t39 * t31 - t38 * t66;
t17 = pkin(3) + t20;
t22 = t38 * t31 + t39 * t66;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t8 = t44 * t17 - t41 * t22;
t6 = -pkin(4) - t8;
t1 = t6 * t26;
t36 = t40 ^ 2;
t62 = mrSges(6,3) * t36;
t9 = t41 * t17 + t44 * t22;
t7 = pkin(8) + t9;
t2 = t7 * t62;
t61 = mrSges(6,3) * t37;
t3 = t7 * t61;
t4 = t8 * mrSges(5,1);
t68 = t1 + t2 + t3 + t4;
t30 = t39 * pkin(2) + pkin(3);
t64 = t38 * pkin(2);
t21 = t44 * t30 - t41 * t64;
t18 = -pkin(4) - t21;
t10 = t18 * t26;
t23 = t41 * t30 + t44 * t64;
t19 = pkin(8) + t23;
t11 = t19 * t62;
t12 = t19 * t61;
t15 = t21 * mrSges(5,1);
t67 = t10 + t11 + t12 + t15;
t65 = pkin(4) * t26;
t63 = t9 * mrSges(5,2);
t60 = t20 * mrSges(4,1);
t59 = t22 * mrSges(4,2);
t58 = t23 * mrSges(5,2);
t57 = Ifges(6,5) * t40 + Ifges(6,6) * t43;
t56 = t36 + t37;
t55 = Ifges(6,2) * t37 + Ifges(5,3) + (Ifges(6,1) * t40 + 0.2e1 * Ifges(6,4) * t43) * t40;
t54 = t56 * t19;
t53 = Ifges(3,3) + Ifges(4,3) + t55;
t52 = t39 * mrSges(4,1) - t38 * mrSges(4,2);
t51 = -mrSges(6,1) * t40 - mrSges(6,2) * t43;
t50 = (t45 * mrSges(3,1) - t42 * mrSges(3,2)) * pkin(1);
t32 = pkin(8) * t62;
t33 = pkin(8) * t61;
t49 = t32 + t33 + t55 - t65;
t5 = [0.2e1 * t60 - 0.2e1 * t59 - 0.2e1 * t63 + Ifges(2,3) + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 + 0.2e1 * t4 + 0.2e1 * t50 + m(6) * (t56 * t7 ^ 2 + t6 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t20 ^ 2 + t22 ^ 2) + m(3) * (t42 ^ 2 + t45 ^ 2) * pkin(1) ^ 2 + t53; t53 + (m(4) * (t20 * t39 + t22 * t38) + t52) * pkin(2) + t50 + m(6) * (t18 * t6 + t7 * t54) + m(5) * (t21 * t8 + t23 * t9) + (-t9 - t23) * mrSges(5,2) + t60 - t59 + t67 + t68; -0.2e1 * t58 + 0.2e1 * t10 + 0.2e1 * t11 + 0.2e1 * t12 + 0.2e1 * t15 + m(6) * (t56 * t19 ^ 2 + t18 ^ 2) + m(5) * (t21 ^ 2 + t23 ^ 2) + t53 + (0.2e1 * t52 + m(4) * (t38 ^ 2 + t39 ^ 2) * pkin(2)) * pkin(2); 0; 0; m(6) * t56 + m(4) + m(5); m(6) * (t56 * t7 * pkin(8) - pkin(4) * t6) - t63 + t49 + t68; m(6) * (-pkin(4) * t18 + pkin(8) * t54) - t58 + t49 + t67; 0; 0.2e1 * t33 + 0.2e1 * t32 - 0.2e1 * t65 + m(6) * (t56 * pkin(8) ^ 2 + pkin(4) ^ 2) + t55; t51 * t7 + t57; t51 * t19 + t57; -t26; t51 * pkin(8) + t57; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
