% Calculate joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:22
% EndTime: 2019-12-31 17:09:23
% DurationCPUTime: 0.38s
% Computational Cost: add. (356->128), mult. (753->194), div. (0->0), fcn. (672->6), ass. (0->60)
t72 = 2 * pkin(5);
t51 = sin(pkin(7));
t52 = cos(pkin(7));
t37 = -t52 * mrSges(4,1) + t51 * mrSges(4,2);
t71 = -m(4) * pkin(2) + t37;
t69 = -t51 / 0.2e1;
t68 = t52 / 0.2e1;
t56 = cos(qJ(2));
t67 = pkin(5) * t56;
t66 = Ifges(4,4) * t51;
t65 = Ifges(4,4) * t52;
t54 = sin(qJ(2));
t64 = t51 * t54;
t63 = t52 * t54;
t62 = pkin(6) + qJ(3);
t53 = sin(qJ(4));
t55 = cos(qJ(4));
t31 = t55 * t51 + t53 * t52;
t22 = t31 * t54;
t30 = -t53 * t51 + t55 * t52;
t23 = t30 * t54;
t61 = -Ifges(5,5) * t23 + Ifges(5,6) * t22;
t24 = mrSges(4,1) * t64 + mrSges(4,2) * t63;
t35 = -t56 * pkin(2) - t54 * qJ(3) - pkin(1);
t16 = t51 * t35 + t52 * t67;
t60 = t51 ^ 2 + t52 ^ 2;
t5 = t22 * mrSges(5,1) + t23 * mrSges(5,2);
t7 = -t30 * mrSges(5,1) + t31 * mrSges(5,2);
t29 = t52 * t35;
t15 = -t51 * t67 + t29;
t59 = -t15 * t51 + t16 * t52;
t58 = pkin(5) ^ 2;
t50 = t56 ^ 2;
t49 = t54 ^ 2;
t46 = t49 * t58;
t44 = -t52 * pkin(3) - pkin(2);
t40 = Ifges(4,1) * t51 + t65;
t39 = Ifges(4,2) * t52 + t66;
t38 = t62 * t52;
t36 = t62 * t51;
t34 = (pkin(3) * t51 + pkin(5)) * t54;
t33 = -t56 * mrSges(4,1) - mrSges(4,3) * t63;
t32 = t56 * mrSges(4,2) - mrSges(4,3) * t64;
t27 = Ifges(5,5) * t31;
t26 = Ifges(5,6) * t30;
t21 = -Ifges(4,5) * t56 + (Ifges(4,1) * t52 - t66) * t54;
t20 = -Ifges(4,6) * t56 + (-Ifges(4,2) * t51 + t65) * t54;
t14 = -t56 * mrSges(5,1) - t23 * mrSges(5,3);
t13 = t56 * mrSges(5,2) - t22 * mrSges(5,3);
t12 = -t53 * t36 + t55 * t38;
t11 = -t55 * t36 - t53 * t38;
t10 = -pkin(6) * t64 + t16;
t9 = Ifges(5,1) * t31 + Ifges(5,4) * t30;
t8 = Ifges(5,4) * t31 + Ifges(5,2) * t30;
t6 = -pkin(6) * t63 + t29 + (-pkin(5) * t51 - pkin(3)) * t56;
t4 = Ifges(5,1) * t23 - Ifges(5,4) * t22 - Ifges(5,5) * t56;
t3 = Ifges(5,4) * t23 - Ifges(5,2) * t22 - Ifges(5,6) * t56;
t2 = t55 * t10 + t53 * t6;
t1 = -t53 * t10 + t55 * t6;
t17 = [0.2e1 * t1 * t14 + 0.2e1 * t2 * t13 + 0.2e1 * t15 * t33 + 0.2e1 * t16 * t32 - t22 * t3 + t23 * t4 + 0.2e1 * t34 * t5 + Ifges(2,3) + (t49 + t50) * mrSges(3,3) * t72 + m(5) * (t1 ^ 2 + t2 ^ 2 + t34 ^ 2) + m(4) * (t15 ^ 2 + t16 ^ 2 + t46) + m(3) * (pkin(1) ^ 2 + t50 * t58 + t46) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(5,3) + Ifges(3,2) + Ifges(4,3)) * t56 + t61) * t56 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t54 + t24 * t72 - t51 * t20 + t52 * t21 + (-Ifges(4,5) * t52 + Ifges(4,6) * t51 + (2 * Ifges(3,4))) * t56) * t54; t20 * t68 + t44 * t5 + t51 * t21 / 0.2e1 + m(5) * (t11 * t1 + t12 * t2 + t44 * t34) + t12 * t13 + t11 * t14 - t22 * t8 / 0.2e1 + t23 * t9 / 0.2e1 - pkin(2) * t24 + t30 * t3 / 0.2e1 + t31 * t4 / 0.2e1 + t34 * t7 + (Ifges(4,5) * t69 - Ifges(4,6) * t52 / 0.2e1 + Ifges(3,6) - t27 / 0.2e1 - t26 / 0.2e1 - pkin(5) * mrSges(3,2)) * t56 + (-t1 * t31 + t2 * t30) * mrSges(5,3) + t59 * mrSges(4,3) + (m(4) * t59 + t52 * t32 - t51 * t33) * qJ(3) + (Ifges(3,5) + t39 * t69 + t40 * t68 + (-mrSges(3,1) + t71) * pkin(5)) * t54; -0.2e1 * pkin(2) * t37 + t30 * t8 + t31 * t9 + t52 * t39 + t51 * t40 + 0.2e1 * t44 * t7 + Ifges(3,3) + m(5) * (t11 ^ 2 + t12 ^ 2 + t44 ^ 2) + m(4) * (t60 * qJ(3) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (-t11 * t31 + t12 * t30) * mrSges(5,3) + 0.2e1 * t60 * qJ(3) * mrSges(4,3); m(4) * t54 * pkin(5) + m(5) * t34 + t24 + t5; m(5) * t44 + t7 + t71; m(4) + m(5); t1 * mrSges(5,1) - t2 * mrSges(5,2) - Ifges(5,3) * t56 - t61; t11 * mrSges(5,1) - t12 * mrSges(5,2) + t26 + t27; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1), t17(2), t17(4), t17(7); t17(2), t17(3), t17(5), t17(8); t17(4), t17(5), t17(6), t17(9); t17(7), t17(8), t17(9), t17(10);];
Mq = res;
