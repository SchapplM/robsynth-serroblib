% Calculate joint inertia matrix for
% S5RPPRR4
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:28
% EndTime: 2020-01-03 11:30:31
% DurationCPUTime: 0.52s
% Computational Cost: add. (683->127), mult. (1427->192), div. (0->0), fcn. (1454->8), ass. (0->58)
t52 = cos(pkin(9));
t72 = t52 ^ 2;
t71 = 2 * qJ(2);
t50 = sin(pkin(9));
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t36 = t57 * t50 + t55 * t52;
t51 = sin(pkin(8));
t31 = t36 * t51;
t35 = -t55 * t50 + t57 * t52;
t32 = t35 * t51;
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t13 = -t56 * t31 - t54 * t32;
t14 = -t54 * t31 + t56 * t32;
t67 = -Ifges(6,5) * t14 - Ifges(6,6) * t13;
t30 = Ifges(5,5) * t32;
t29 = Ifges(5,6) * t31;
t70 = t50 * t51;
t69 = t52 * t51;
t68 = -Ifges(6,3) - Ifges(5,3);
t53 = cos(pkin(8));
t40 = -t53 * pkin(2) - t51 * qJ(3) - pkin(1);
t34 = t52 * t40;
t18 = -pkin(6) * t69 + t34 + (-qJ(2) * t50 - pkin(3)) * t53;
t65 = qJ(2) * t53;
t27 = t50 * t40 + t52 * t65;
t23 = -pkin(6) * t70 + t27;
t7 = t55 * t18 + t57 * t23;
t66 = mrSges(4,1) * t70 + mrSges(4,2) * t69;
t46 = t51 * qJ(2);
t39 = pkin(3) * t70 + t46;
t64 = t31 * mrSges(5,1) + t32 * mrSges(5,2);
t63 = -t13 * mrSges(6,1) + t14 * mrSges(6,2);
t21 = t56 * t35 - t54 * t36;
t22 = t54 * t35 + t56 * t36;
t62 = t21 * mrSges(6,1) - t22 * mrSges(6,2);
t6 = t57 * t18 - t55 * t23;
t4 = -t53 * pkin(4) - t32 * pkin(7) + t6;
t5 = -t31 * pkin(7) + t7;
t2 = t56 * t4 - t54 * t5;
t3 = t54 * t4 + t56 * t5;
t61 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t67;
t60 = (t56 * mrSges(6,1) - t54 * mrSges(6,2)) * pkin(4);
t58 = qJ(2) ^ 2;
t49 = t53 ^ 2;
t48 = t51 ^ 2;
t47 = t51 * mrSges(3,2);
t45 = t48 * t58;
t38 = -t53 * mrSges(4,1) - mrSges(4,3) * t69;
t37 = t53 * mrSges(4,2) - mrSges(4,3) * t70;
t26 = -t50 * t65 + t34;
t25 = -t53 * mrSges(5,1) - t32 * mrSges(5,3);
t24 = t53 * mrSges(5,2) - t31 * mrSges(5,3);
t20 = t31 * pkin(4) + t39;
t9 = -t53 * mrSges(6,1) - t14 * mrSges(6,3);
t8 = t53 * mrSges(6,2) + t13 * mrSges(6,3);
t1 = [-0.2e1 * pkin(1) * t47 + Ifges(5,2) * t31 ^ 2 + Ifges(6,2) * t13 ^ 2 + 0.2e1 * t26 * t38 + 0.2e1 * t39 * t64 + 0.2e1 * t27 * t37 + 0.2e1 * t20 * t63 + 0.2e1 * t7 * t24 + 0.2e1 * t6 * t25 + 0.2e1 * t3 * t8 + 0.2e1 * t2 * t9 + Ifges(2,3) + (Ifges(5,1) * t32 - 0.2e1 * Ifges(5,4) * t31) * t32 + (Ifges(6,1) * t14 + 0.2e1 * Ifges(6,4) * t13) * t14 + (t48 + t49) * mrSges(3,3) * t71 + m(6) * (t2 ^ 2 + t20 ^ 2 + t3 ^ 2) + m(5) * (t39 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2 + t45) + m(3) * (pkin(1) ^ 2 + t49 * t58 + t45) + (0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * t30 + 0.2e1 * t29 + (Ifges(4,3) + Ifges(3,2) - t68) * t53 + 0.2e1 * t67) * t53 + (t66 * t71 + (Ifges(4,1) * t72 + Ifges(3,1) + (-0.2e1 * Ifges(4,4) * t52 + Ifges(4,2) * t50) * t50) * t51 + 0.2e1 * (-Ifges(4,5) * t52 + Ifges(4,6) * t50 + Ifges(3,4)) * t53) * t51; -m(3) * pkin(1) - t53 * mrSges(3,1) + t21 * t9 + t22 * t8 + t36 * t24 + t35 * t25 + t50 * t37 + t52 * t38 + t47 + m(6) * (t21 * t2 + t22 * t3) + m(5) * (t35 * t6 + t36 * t7) + m(4) * (t52 * t26 + t50 * t27); m(3) + m(4) * (t50 ^ 2 + t72) + m(5) * (t35 ^ 2 + t36 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2); m(4) * t46 + m(5) * t39 + m(6) * t20 + t63 + t64 + t66; 0; m(4) + m(5) + m(6); t6 * mrSges(5,1) - t7 * mrSges(5,2) - t29 + t30 + t68 * t53 + (m(6) * (t2 * t56 + t3 * t54) + t54 * t8 + t56 * t9) * pkin(4) + t61; t35 * mrSges(5,1) - t36 * mrSges(5,2) + m(6) * (t21 * t56 + t22 * t54) * pkin(4) + t62; 0; m(6) * (t54 ^ 2 + t56 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t60 - t68; -Ifges(6,3) * t53 + t61; t62; 0; Ifges(6,3) + t60; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
