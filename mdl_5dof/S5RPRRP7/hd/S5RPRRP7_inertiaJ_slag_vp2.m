% Calculate joint inertia matrix for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:34
% EndTime: 2019-12-31 18:44:36
% DurationCPUTime: 0.50s
% Computational Cost: add. (395->154), mult. (799->208), div. (0->0), fcn. (573->6), ass. (0->67)
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t65 = t49 ^ 2 + t51 ^ 2;
t47 = sin(pkin(8));
t35 = pkin(1) * t47 + pkin(6);
t83 = 0.2e1 * t35;
t82 = (mrSges(6,2) + mrSges(5,3)) * t65;
t50 = sin(qJ(3));
t81 = t50 * pkin(7);
t52 = cos(qJ(3));
t80 = t52 * pkin(3);
t79 = Ifges(5,4) * t49;
t78 = Ifges(5,4) * t51;
t77 = Ifges(6,5) * t49;
t76 = Ifges(6,5) * t51;
t75 = Ifges(5,6) * t52;
t74 = Ifges(6,6) * t52;
t73 = t35 * t52;
t72 = t49 * t50;
t71 = t50 * t51;
t48 = cos(pkin(8));
t36 = -pkin(1) * t48 - pkin(2);
t14 = t36 - t80 - t81;
t70 = t51 * t14;
t22 = -mrSges(5,1) * t51 + mrSges(5,2) * t49;
t69 = mrSges(4,1) - t22;
t68 = Ifges(6,2) + Ifges(5,3);
t4 = t49 * t14 + t51 * t73;
t17 = t52 * mrSges(6,1) + mrSges(6,2) * t71;
t67 = t65 * t81;
t66 = t65 * pkin(7) ^ 2;
t44 = t50 ^ 2;
t46 = t52 ^ 2;
t64 = t44 + t46;
t63 = -Ifges(6,6) * t72 + (-Ifges(6,4) - Ifges(5,5)) * t71;
t58 = t49 * mrSges(6,1) - t51 * mrSges(6,3);
t12 = t58 * t50;
t57 = -pkin(4) * t49 + qJ(5) * t51;
t5 = (t35 - t57) * t50;
t62 = m(6) * t5 + t12;
t3 = -t49 * t73 + t70;
t60 = -t3 * t49 + t4 * t51;
t59 = t49 * mrSges(5,1) + t51 * mrSges(5,2);
t1 = -qJ(5) * t52 + t4;
t15 = mrSges(5,2) * t52 - mrSges(5,3) * t72;
t16 = -mrSges(5,1) * t52 - mrSges(5,3) * t71;
t18 = -mrSges(6,2) * t72 - mrSges(6,3) * t52;
t2 = -t70 + (t35 * t49 + pkin(4)) * t52;
t56 = m(6) * (t1 * t51 + t2 * t49) + (t15 + t18) * t51 + (-t16 + t17) * t49;
t55 = m(6) * t57 - t58 - t59;
t40 = Ifges(6,4) * t49;
t39 = Ifges(5,5) * t49;
t38 = Ifges(5,6) * t51;
t34 = t35 ^ 2;
t27 = t44 * t34;
t26 = Ifges(5,1) * t49 + t78;
t25 = Ifges(6,1) * t49 - t76;
t24 = Ifges(5,2) * t51 + t79;
t23 = -Ifges(6,3) * t51 + t77;
t21 = -mrSges(6,1) * t51 - mrSges(6,3) * t49;
t20 = -pkin(4) * t51 - qJ(5) * t49 - pkin(3);
t13 = t59 * t50;
t9 = -Ifges(5,5) * t52 + (Ifges(5,1) * t51 - t79) * t50;
t8 = -Ifges(6,4) * t52 + (Ifges(6,1) * t51 + t77) * t50;
t7 = -t75 + (-Ifges(5,2) * t49 + t78) * t50;
t6 = -t74 + (Ifges(6,3) * t49 + t76) * t50;
t10 = [0.2e1 * t1 * t18 + 0.2e1 * t5 * t12 + 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + Ifges(2,3) + Ifges(3,3) + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t27) + m(4) * (t34 * t46 + t36 ^ 2 + t27) + m(3) * (t47 ^ 2 + t48 ^ 2) * pkin(1) ^ 2 + (-0.2e1 * t36 * mrSges(4,1) + (Ifges(4,2) + t68) * t52 + t63) * t52 + (0.2e1 * t36 * mrSges(4,2) + Ifges(4,1) * t50 + 0.2e1 * Ifges(4,4) * t52 + t13 * t83 + (t8 + t9) * t51 + (t6 - t7 + t75) * t49) * t50 + 0.2e1 * (t48 * mrSges(3,1) - t47 * mrSges(3,2)) * pkin(1) + t64 * mrSges(4,3) * t83; (-t13 - t62) * t52 + (m(5) * (t60 - t73) + t56) * t50; m(3) + m(4) * t64 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t44 * t65 + t46); -pkin(3) * t13 + t5 * t21 + t62 * t20 + (-t40 / 0.2e1 - t39 / 0.2e1 - t38 / 0.2e1 + Ifges(4,6) - t35 * mrSges(4,2)) * t52 + (-t6 / 0.2e1 + t7 / 0.2e1 + t74 / 0.2e1 + t1 * mrSges(6,2) + t4 * mrSges(5,3)) * t51 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t3 * mrSges(5,3)) * t49 + (m(5) * t60 + t56) * pkin(7) + (Ifges(4,5) + (t25 / 0.2e1 + t26 / 0.2e1) * t51 + (t23 / 0.2e1 - t24 / 0.2e1) * t49 + (-m(5) * pkin(3) - t69) * t35) * t50; (-t21 + t69) * t52 + m(5) * (t67 + t80) + m(6) * (-t20 * t52 + t67) + (-mrSges(4,2) + t82) * t50; -0.2e1 * pkin(3) * t22 + 0.2e1 * t20 * t21 + Ifges(4,3) + (-t23 + t24) * t51 + (t25 + t26) * t49 + m(6) * (t20 ^ 2 + t66) + m(5) * (pkin(3) ^ 2 + t66) + 0.2e1 * pkin(7) * t82; -Ifges(5,6) * t72 + m(6) * (-pkin(4) * t2 + qJ(5) * t1) - t4 * mrSges(5,2) + t3 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t17 + t1 * mrSges(6,3) + qJ(5) * t18 - t68 * t52 - t63; t55 * t50; mrSges(6,2) * t57 - Ifges(6,6) * t51 + pkin(7) * t55 + t38 + t39 + t40; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t68; m(6) * t2 + t17; m(6) * t72; (m(6) * pkin(7) + mrSges(6,2)) * t49; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
