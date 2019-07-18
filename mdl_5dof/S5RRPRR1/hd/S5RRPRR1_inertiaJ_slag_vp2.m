% Calculate joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:45
% EndTime: 2019-07-18 17:20:47
% DurationCPUTime: 0.39s
% Computational Cost: add. (449->122), mult. (881->169), div. (0->0), fcn. (782->6), ass. (0->59)
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t30 = -t50 * mrSges(6,1) + t47 * mrSges(6,2);
t83 = mrSges(5,1) - t30;
t49 = sin(qJ(2));
t43 = t49 ^ 2;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t53 = pkin(2) + pkin(1);
t82 = (-t48 * mrSges(5,2) + t83 * t51) * t53;
t52 = cos(qJ(2));
t70 = pkin(3) + qJ(3);
t31 = t70 * t52;
t62 = t70 * t49;
t15 = t48 * t31 + t51 * t62;
t81 = t15 ^ 2;
t80 = 0.2e1 * t15;
t34 = t53 * t52;
t79 = -0.2e1 * t34;
t17 = t51 * t31 - t48 * t62;
t29 = t48 * t52 + t51 * t49;
t19 = -t29 * pkin(4) - t34;
t3 = t50 * t17 + t47 * t19;
t78 = t3 * t50;
t77 = Ifges(6,4) * t47;
t76 = Ifges(6,4) * t50;
t75 = t15 * t51;
t74 = t29 * t47;
t73 = t29 * t50;
t72 = t47 * mrSges(6,3);
t28 = t48 * t49 - t51 * t52;
t69 = Ifges(6,5) * t73 + Ifges(6,3) * t28;
t68 = Ifges(6,5) * t47 + Ifges(6,6) * t50;
t67 = t47 ^ 2 + t50 ^ 2;
t66 = 2 * pkin(1) * mrSges(4,1);
t32 = Ifges(6,2) * t50 + t77;
t33 = Ifges(6,1) * t47 + t76;
t65 = t50 * t32 + t47 * t33 + Ifges(5,3);
t64 = m(6) * t67;
t36 = t48 * t53 + pkin(4);
t61 = t67 * t36;
t60 = mrSges(6,1) * t47 + mrSges(6,2) * t50;
t59 = 0.2e1 * t67 * mrSges(6,3);
t10 = -t28 * mrSges(6,2) - t29 * t72;
t11 = t28 * mrSges(6,1) - mrSges(6,3) * t73;
t2 = -t47 * t17 + t50 * t19;
t58 = m(6) * (-t2 * t47 + t78) + t50 * t10 - t47 * t11;
t6 = Ifges(6,6) * t28 + (-Ifges(6,2) * t47 + t76) * t29;
t7 = Ifges(6,5) * t28 + (Ifges(6,1) * t50 - t77) * t29;
t57 = -t17 * mrSges(5,2) + mrSges(6,3) * t78 - t2 * t72 - t32 * t74 / 0.2e1 + t33 * t73 / 0.2e1 + Ifges(5,5) * t29 + t47 * t7 / 0.2e1 + t50 * t6 / 0.2e1 + (t68 / 0.2e1 - Ifges(5,6)) * t28 - t83 * t15;
t56 = pkin(1) ^ 2;
t54 = qJ(3) ^ 2;
t46 = t53 ^ 2;
t45 = t52 ^ 2;
t39 = t49 * mrSges(4,2);
t38 = t51 ^ 2 * t46;
t21 = t29 * mrSges(5,2);
t9 = t60 * t29;
t1 = [0.2e1 * t3 * t10 + 0.2e1 * t2 * t11 + t9 * t80 + t21 * t79 + Ifges(2,3) + (mrSges(5,1) * t79 - 0.2e1 * t17 * mrSges(5,3) + Ifges(5,2) * t28 + t69) * t28 + m(6) * (t2 ^ 2 + t3 ^ 2 + t81) + m(5) * (t17 ^ 2 + t34 ^ 2 + t81) + m(4) * (t43 * t54 + (t54 + t56) * t45) + (Ifges(4,1) + Ifges(3,1)) * t43 + (mrSges(5,3) * t80 + Ifges(5,1) * t29 - t47 * t6 + t50 * t7 + (-Ifges(6,6) * t47 - (2 * Ifges(5,4))) * t28) * t29 + 0.2e1 * (t43 + t45) * qJ(3) * mrSges(4,3) + (-0.2e1 * pkin(1) * t39 + 0.2e1 * (Ifges(3,4) + Ifges(4,4)) * t49 + (Ifges(3,2) + Ifges(4,2) + t66) * t52) * t52; t57 + t58 * t36 + (-t51 * t9 + (-t48 * t28 - t51 * t29) * mrSges(5,3) + m(5) * (t17 * t48 - t75) - m(6) * t75) * t53 + (-qJ(3) * mrSges(4,1) + Ifges(3,5) + Ifges(4,5) + (-m(4) * qJ(3) - mrSges(4,3)) * pkin(1)) * t49 + (-qJ(3) * mrSges(4,2) + Ifges(3,6) + Ifges(4,6)) * t52; m(4) * t56 + t66 + Ifges(3,3) + Ifges(4,3) + t36 * t59 + m(6) * (t67 * t36 ^ 2 + t38) + m(5) * (t48 ^ 2 * t46 + t38) + 0.2e1 * t82 + t65; t28 * mrSges(5,1) + t47 * t10 + t50 * t11 + t21 + t39 + (-m(4) * pkin(1) - mrSges(4,1)) * t52 + m(6) * (t50 * t2 + t47 * t3) - m(5) * t34; 0; m(4) + m(5) + t64; t58 * pkin(4) + t57; m(6) * pkin(4) * t61 + t82 + (t67 * pkin(4) + t61) * mrSges(6,3) + t65; 0; t65 + (t64 * pkin(4) + t59) * pkin(4); t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t74 + t69; -t60 * t36 + t68; -t30; -t60 * pkin(4) + t68; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
