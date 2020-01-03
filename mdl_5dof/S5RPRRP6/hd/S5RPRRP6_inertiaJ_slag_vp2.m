% Calculate joint inertia matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:10
% EndTime: 2019-12-31 18:42:12
% DurationCPUTime: 0.61s
% Computational Cost: add. (389->159), mult. (779->216), div. (0->0), fcn. (577->6), ass. (0->65)
t84 = Ifges(5,5) + Ifges(6,5);
t65 = Ifges(5,6) + Ifges(6,6);
t83 = -2 * mrSges(6,3);
t47 = sin(pkin(8));
t35 = t47 * pkin(1) + pkin(6);
t82 = 0.2e1 * t35;
t49 = sin(qJ(4));
t50 = sin(qJ(3));
t11 = (pkin(4) * t49 + t35) * t50;
t51 = cos(qJ(4));
t67 = t50 * t51;
t68 = t49 * t50;
t13 = mrSges(6,1) * t68 + mrSges(6,2) * t67;
t81 = m(6) * t11 + t13;
t80 = t84 * t49 + t65 * t51;
t79 = m(6) * pkin(4);
t76 = t50 * pkin(7);
t52 = cos(qJ(3));
t75 = t52 * pkin(3);
t48 = cos(pkin(8));
t36 = -t48 * pkin(1) - pkin(2);
t15 = t36 - t75 - t76;
t69 = t35 * t52;
t4 = t49 * t15 + t51 * t69;
t74 = mrSges(5,2) * t51;
t73 = Ifges(5,4) * t49;
t72 = Ifges(5,4) * t51;
t71 = Ifges(6,4) * t49;
t70 = Ifges(6,4) * t51;
t23 = -t51 * mrSges(5,1) + t49 * mrSges(5,2);
t66 = mrSges(4,1) - t23;
t64 = Ifges(5,3) + Ifges(6,3);
t63 = -qJ(5) - pkin(7);
t62 = t84 * t67;
t59 = t49 ^ 2 + t51 ^ 2;
t44 = t50 ^ 2;
t46 = t52 ^ 2;
t58 = t44 + t46;
t57 = qJ(5) * t50;
t56 = t59 * mrSges(5,3);
t22 = -t51 * mrSges(6,1) + t49 * mrSges(6,2);
t55 = mrSges(5,1) * t49 + t74;
t37 = -t51 * pkin(4) - pkin(3);
t34 = t35 ^ 2;
t29 = t44 * t34;
t28 = Ifges(5,1) * t49 + t72;
t27 = Ifges(6,1) * t49 + t70;
t26 = Ifges(5,2) * t51 + t73;
t25 = Ifges(6,2) * t51 + t71;
t24 = t63 * t51;
t21 = t63 * t49;
t19 = -t52 * mrSges(5,1) - mrSges(5,3) * t67;
t18 = -t52 * mrSges(6,1) - mrSges(6,3) * t67;
t17 = t52 * mrSges(5,2) - mrSges(5,3) * t68;
t16 = t52 * mrSges(6,2) - mrSges(6,3) * t68;
t14 = t55 * t50;
t10 = t51 * t15;
t8 = -Ifges(5,5) * t52 + (Ifges(5,1) * t51 - t73) * t50;
t7 = -Ifges(6,5) * t52 + (Ifges(6,1) * t51 - t71) * t50;
t6 = -Ifges(5,6) * t52 + (-Ifges(5,2) * t49 + t72) * t50;
t5 = -Ifges(6,6) * t52 + (-Ifges(6,2) * t49 + t70) * t50;
t3 = -t49 * t69 + t10;
t2 = -t49 * t57 + t4;
t1 = -t51 * t57 + t10 + (-t35 * t49 - pkin(4)) * t52;
t9 = [0.2e1 * t1 * t18 + 0.2e1 * t11 * t13 + 0.2e1 * t2 * t16 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t19 + Ifges(2,3) + Ifges(3,3) + m(6) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t29) + m(4) * (t46 * t34 + t36 ^ 2 + t29) + m(3) * (t47 ^ 2 + t48 ^ 2) * pkin(1) ^ 2 + (-0.2e1 * t36 * mrSges(4,1) + (Ifges(4,2) + t64) * t52 - t62) * t52 + (0.2e1 * t36 * mrSges(4,2) + Ifges(4,1) * t50 + 0.2e1 * Ifges(4,4) * t52 + t14 * t82 + (t7 + t8) * t51 + (t52 * t65 - t5 - t6) * t49) * t50 + 0.2e1 * (t48 * mrSges(3,1) - t47 * mrSges(3,2)) * pkin(1) + t58 * mrSges(4,3) * t82; (-t14 - t81) * t52 + ((t16 + t17) * t51 + (-t18 - t19) * t49 + m(5) * (-t3 * t49 + t4 * t51 - t69) + m(6) * (-t1 * t49 + t2 * t51)) * t50; m(3) + m(4) * t58 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t44 * t59 + t46); t37 * t13 + t21 * t18 + t11 * t22 - t24 * t16 - pkin(3) * t14 + m(6) * (t21 * t1 + t37 * t11 - t24 * t2) + (t5 / 0.2e1 + t6 / 0.2e1 + t2 * mrSges(6,3) + t4 * mrSges(5,3) + (m(5) * t4 + t17) * pkin(7)) * t51 + (t7 / 0.2e1 + t8 / 0.2e1 - t1 * mrSges(6,3) - t3 * mrSges(5,3) + (-m(5) * t3 - t19) * pkin(7)) * t49 + (Ifges(4,5) + (-m(5) * pkin(3) - t66) * t35 + (t27 / 0.2e1 + t28 / 0.2e1) * t51 + (-t25 / 0.2e1 - t26 / 0.2e1) * t49) * t50 + (Ifges(4,6) - t35 * mrSges(4,2) - t80 / 0.2e1) * t52; (-t22 + t66) * t52 + (t59 * mrSges(6,3) - mrSges(4,2) + t56) * t50 + m(5) * (t59 * t76 + t75) + m(6) * (-t37 * t52 + (-t21 * t49 - t24 * t51) * t50); -0.2e1 * pkin(3) * t23 + 0.2e1 * t37 * t22 + Ifges(4,3) + 0.2e1 * pkin(7) * t56 + m(5) * (pkin(7) ^ 2 * t59 + pkin(3) ^ 2) + m(6) * (t21 ^ 2 + t24 ^ 2 + t37 ^ 2) + (t24 * t83 + t25 + t26) * t51 + (t21 * t83 + t27 + t28) * t49; t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) - t64 * t52 - t65 * t68 + (m(6) * t1 + t18) * pkin(4) + t62; (-t74 + (-mrSges(5,1) - t79) * t49) * t50 - t13; t21 * mrSges(6,1) + t24 * mrSges(6,2) - t55 * pkin(7) + (m(6) * t21 - t49 * mrSges(6,3)) * pkin(4) + t80; (0.2e1 * mrSges(6,1) + t79) * pkin(4) + t64; t81; -m(6) * t52; m(6) * t37 + t22; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
