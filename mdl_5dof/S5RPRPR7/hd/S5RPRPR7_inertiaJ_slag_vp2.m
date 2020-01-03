% Calculate joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:55
% EndTime: 2019-12-31 18:18:56
% DurationCPUTime: 0.39s
% Computational Cost: add. (506->127), mult. (955->194), div. (0->0), fcn. (897->8), ass. (0->54)
t45 = cos(qJ(3));
t71 = t45 ^ 2;
t38 = sin(pkin(9));
t40 = cos(pkin(9));
t43 = sin(qJ(3));
t18 = t38 * t43 - t40 * t45;
t20 = t38 * t45 + t40 * t43;
t42 = sin(qJ(5));
t61 = t20 * t42;
t10 = -t18 * mrSges(6,2) - mrSges(6,3) * t61;
t44 = cos(qJ(5));
t60 = t20 * t44;
t11 = t18 * mrSges(6,1) - mrSges(6,3) * t60;
t70 = t44 * t10 - t42 * t11;
t39 = sin(pkin(8));
t29 = t39 * pkin(1) + pkin(6);
t53 = qJ(4) + t29;
t16 = t53 * t45;
t51 = t53 * t43;
t6 = t38 * t16 + t40 * t51;
t69 = t6 ^ 2;
t68 = 0.2e1 * t6;
t67 = t18 ^ 2;
t41 = cos(pkin(8));
t31 = -t41 * pkin(1) - pkin(2);
t21 = -t45 * pkin(3) + t31;
t66 = 0.2e1 * t21;
t65 = t44 / 0.2e1;
t64 = t6 * t18;
t63 = Ifges(6,4) * t42;
t62 = Ifges(6,4) * t44;
t57 = Ifges(6,5) * t60 + Ifges(6,3) * t18;
t56 = Ifges(6,5) * t42 + Ifges(6,6) * t44;
t55 = t42 ^ 2 + t44 ^ 2;
t54 = t43 ^ 2 + t71;
t28 = t38 * pkin(3) + pkin(7);
t52 = t55 * t28;
t5 = t18 * pkin(4) - t20 * pkin(7) + t21;
t8 = t40 * t16 - t38 * t51;
t1 = -t42 * t8 + t44 * t5;
t2 = t42 * t5 + t44 * t8;
t50 = -t1 * t42 + t2 * t44;
t49 = -t45 * mrSges(4,1) + t43 * mrSges(4,2);
t48 = mrSges(6,1) * t42 + mrSges(6,2) * t44;
t30 = -t40 * pkin(3) - pkin(4);
t24 = Ifges(6,1) * t42 + t62;
t23 = Ifges(6,2) * t44 + t63;
t22 = -t44 * mrSges(6,1) + t42 * mrSges(6,2);
t17 = t20 ^ 2;
t14 = t20 * mrSges(5,2);
t9 = t48 * t20;
t4 = Ifges(6,5) * t18 + (Ifges(6,1) * t44 - t63) * t20;
t3 = Ifges(6,6) * t18 + (-Ifges(6,2) * t42 + t62) * t20;
t7 = [Ifges(2,3) + Ifges(3,3) + t9 * t68 + 0.2e1 * t1 * t11 + 0.2e1 * t2 * t10 + Ifges(4,2) * t71 + 0.2e1 * t31 * t49 + t14 * t66 + (mrSges(5,1) * t66 - 0.2e1 * t8 * mrSges(5,3) + Ifges(5,2) * t18 + t57) * t18 + (mrSges(5,3) * t68 + Ifges(5,1) * t20 - t42 * t3 + t44 * t4 + (-Ifges(6,6) * t42 - (2 * Ifges(5,4))) * t18) * t20 + m(6) * (t1 ^ 2 + t2 ^ 2 + t69) + m(5) * (t21 ^ 2 + t8 ^ 2 + t69) + m(4) * (t54 * t29 ^ 2 + t31 ^ 2) + m(3) * (t39 ^ 2 + t41 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t43 + 0.2e1 * Ifges(4,4) * t45) * t43 + 0.2e1 * (t41 * mrSges(3,1) - t39 * mrSges(3,2)) * pkin(1) + 0.2e1 * t54 * t29 * mrSges(4,3); t18 * t9 + t70 * t20 + m(6) * (t50 * t20 + t64) + m(5) * (t8 * t20 + t64); m(3) + m(6) * (t55 * t17 + t67) + m(4) * t54 + m(5) * (t17 + t67); -Ifges(5,6) * t18 + Ifges(4,5) * t43 + Ifges(4,6) * t45 + t42 * t4 / 0.2e1 + t3 * t65 + t6 * t22 + t18 * t56 / 0.2e1 - t8 * mrSges(5,2) - t6 * mrSges(5,1) + (-t43 * mrSges(4,1) - t45 * mrSges(4,2)) * t29 + t50 * mrSges(6,3) + (Ifges(5,5) + t24 * t65 - t42 * t23 / 0.2e1) * t20 + (m(5) * (t38 * t8 - t40 * t6) + (-t38 * t18 - t40 * t20) * mrSges(5,3)) * pkin(3) + (m(6) * t6 + t9) * t30 + (m(6) * t50 + t70) * t28; -t14 + (-mrSges(5,1) + t22) * t18 + t55 * t20 * mrSges(6,3) + m(6) * (t30 * t18 + t20 * t52) + m(5) * (-t18 * t40 + t20 * t38) * pkin(3) - t49; 0.2e1 * t30 * t22 + t44 * t23 + t42 * t24 + Ifges(4,3) + Ifges(5,3) + m(6) * (t55 * t28 ^ 2 + t30 ^ 2) + m(5) * (t38 ^ 2 + t40 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (t40 * mrSges(5,1) - t38 * mrSges(5,2)) * pkin(3) + 0.2e1 * mrSges(6,3) * t52; t18 * mrSges(5,1) + t42 * t10 + t44 * t11 + t14 + m(6) * (t44 * t1 + t42 * t2) + m(5) * t21; 0; 0; m(6) * t55 + m(5); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,6) * t61 + t57; -t9; -t48 * t28 + t56; -t22; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
