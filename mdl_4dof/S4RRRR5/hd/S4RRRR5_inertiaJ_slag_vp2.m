% Calculate joint inertia matrix for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:23
% EndTime: 2019-12-31 17:27:24
% DurationCPUTime: 0.42s
% Computational Cost: add. (438->143), mult. (914->219), div. (0->0), fcn. (809->6), ass. (0->66)
t75 = 2 * pkin(5);
t55 = cos(qJ(3));
t74 = t55 / 0.2e1;
t73 = -pkin(7) - pkin(6);
t56 = cos(qJ(2));
t72 = pkin(5) * t56;
t52 = sin(qJ(3));
t71 = Ifges(4,4) * t52;
t70 = Ifges(4,4) * t55;
t53 = sin(qJ(2));
t69 = t52 * t53;
t68 = t53 * t55;
t67 = -Ifges(5,3) - Ifges(4,3);
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t31 = t51 * t55 + t54 * t52;
t23 = t31 * t53;
t30 = -t51 * t52 + t54 * t55;
t24 = t30 * t53;
t66 = -Ifges(5,5) * t24 + Ifges(5,6) * t23;
t35 = -t56 * pkin(2) - t53 * pkin(6) - pkin(1);
t18 = t52 * t35 + t55 * t72;
t65 = t52 ^ 2 + t55 ^ 2;
t64 = t52 * mrSges(4,1) + t55 * mrSges(4,2);
t29 = t55 * t35;
t17 = -t52 * t72 + t29;
t63 = -t17 * t52 + t18 * t55;
t11 = -pkin(7) * t69 + t18;
t7 = -pkin(7) * t68 + t29 + (-pkin(5) * t52 - pkin(3)) * t56;
t2 = -t51 * t11 + t54 * t7;
t3 = t54 * t11 + t51 * t7;
t62 = t2 * mrSges(5,1) - t3 * mrSges(5,2) - t66;
t39 = t73 * t52;
t40 = t73 * t55;
t13 = t54 * t39 + t51 * t40;
t14 = t51 * t39 - t54 * t40;
t26 = Ifges(5,6) * t30;
t27 = Ifges(5,5) * t31;
t61 = t13 * mrSges(5,1) - t14 * mrSges(5,2) + t26 + t27;
t60 = (t54 * mrSges(5,1) - t51 * mrSges(5,2)) * pkin(3);
t58 = pkin(5) ^ 2;
t50 = t56 ^ 2;
t48 = t53 ^ 2;
t46 = t48 * t58;
t45 = Ifges(4,5) * t52;
t44 = Ifges(4,6) * t55;
t43 = -t55 * pkin(3) - pkin(2);
t41 = Ifges(4,5) * t68;
t38 = Ifges(4,1) * t52 + t70;
t37 = Ifges(4,2) * t55 + t71;
t36 = -t55 * mrSges(4,1) + t52 * mrSges(4,2);
t34 = (pkin(3) * t52 + pkin(5)) * t53;
t33 = -t56 * mrSges(4,1) - mrSges(4,3) * t68;
t32 = t56 * mrSges(4,2) - mrSges(4,3) * t69;
t25 = t64 * t53;
t22 = -Ifges(4,5) * t56 + (Ifges(4,1) * t55 - t71) * t53;
t21 = -Ifges(4,6) * t56 + (-Ifges(4,2) * t52 + t70) * t53;
t16 = -t56 * mrSges(5,1) - t24 * mrSges(5,3);
t15 = t56 * mrSges(5,2) - t23 * mrSges(5,3);
t10 = Ifges(5,1) * t31 + Ifges(5,4) * t30;
t9 = Ifges(5,4) * t31 + Ifges(5,2) * t30;
t8 = -t30 * mrSges(5,1) + t31 * mrSges(5,2);
t6 = t23 * mrSges(5,1) + t24 * mrSges(5,2);
t5 = Ifges(5,1) * t24 - Ifges(5,4) * t23 - Ifges(5,5) * t56;
t4 = Ifges(5,4) * t24 - Ifges(5,2) * t23 - Ifges(5,6) * t56;
t1 = [0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + 0.2e1 * t17 * t33 + 0.2e1 * t18 * t32 - t23 * t4 + t24 * t5 + 0.2e1 * t34 * t6 + Ifges(2,3) + (t48 + t50) * mrSges(3,3) * t75 + m(5) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2 + t46) + m(3) * (pkin(1) ^ 2 + t50 * t58 + t46) + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t53 - t52 * t21 + t55 * t22 + t25 * t75) * t53 + (0.2e1 * pkin(1) * mrSges(3,1) - t41 + (Ifges(3,2) - t67) * t56 + (Ifges(4,6) * t52 + (2 * Ifges(3,4))) * t53 + t66) * t56; t21 * t74 + t52 * t22 / 0.2e1 + t34 * t8 + t43 * t6 - t23 * t9 / 0.2e1 + t24 * t10 / 0.2e1 - pkin(2) * t25 + t30 * t4 / 0.2e1 + t31 * t5 / 0.2e1 + t14 * t15 + t13 * t16 + m(5) * (t13 * t2 + t14 * t3 + t43 * t34) + (-t45 / 0.2e1 - t44 / 0.2e1 - t27 / 0.2e1 - t26 / 0.2e1 + Ifges(3,6) - pkin(5) * mrSges(3,2)) * t56 + (-t2 * t31 + t3 * t30) * mrSges(5,3) + t63 * mrSges(4,3) + (m(4) * t63 + t55 * t32 - t52 * t33) * pkin(6) + (Ifges(3,5) - t52 * t37 / 0.2e1 + t38 * t74 + (-m(4) * pkin(2) - mrSges(3,1) + t36) * pkin(5)) * t53; -0.2e1 * pkin(2) * t36 + t31 * t10 + t30 * t9 + t55 * t37 + t52 * t38 + 0.2e1 * t43 * t8 + Ifges(3,3) + m(5) * (t13 ^ 2 + t14 ^ 2 + t43 ^ 2) + m(4) * (t65 * pkin(6) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (-t13 * t31 + t14 * t30) * mrSges(5,3) + 0.2e1 * t65 * pkin(6) * mrSges(4,3); -Ifges(4,6) * t69 + t17 * mrSges(4,1) - t18 * mrSges(4,2) + t41 + t67 * t56 + (m(5) * (t2 * t54 + t3 * t51) + t51 * t15 + t54 * t16) * pkin(3) + t62; t44 + t45 - t64 * pkin(6) + (m(5) * (t13 * t54 + t14 * t51) + (t51 * t30 - t54 * t31) * mrSges(5,3)) * pkin(3) + t61; m(5) * (t51 ^ 2 + t54 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t60 - t67; -Ifges(5,3) * t56 + t62; t61; Ifges(5,3) + t60; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
