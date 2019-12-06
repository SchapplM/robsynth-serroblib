% Calculate joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:39
% EndTime: 2019-12-05 15:18:41
% DurationCPUTime: 0.39s
% Computational Cost: add. (427->140), mult. (1174->218), div. (0->0), fcn. (1250->12), ass. (0->70)
t82 = 2 * pkin(8);
t81 = m(6) * pkin(9) + mrSges(6,3);
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t25 = -t51 * mrSges(6,1) + t48 * mrSges(6,2);
t80 = -m(6) * pkin(4) - mrSges(5,1) + t25;
t43 = sin(pkin(6));
t44 = sin(pkin(5));
t45 = cos(pkin(11));
t46 = cos(pkin(6));
t47 = cos(pkin(5));
t16 = -t44 * t45 * t43 + t47 * t46;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t42 = sin(pkin(11));
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t66 = t45 * t46;
t68 = t43 * t50;
t8 = t47 * t68 + (t42 * t53 + t50 * t66) * t44;
t3 = -t16 * t52 + t8 * t49;
t79 = t3 ^ 2;
t67 = t43 * t53;
t6 = -t47 * t67 + (t42 * t50 - t53 * t66) * t44;
t78 = t6 ^ 2;
t17 = -t52 * t46 + t49 * t68;
t77 = t17 ^ 2;
t76 = t51 / 0.2e1;
t75 = pkin(8) * t52;
t74 = t17 * t3;
t73 = t3 * t49;
t72 = Ifges(6,4) * t48;
t71 = Ifges(6,4) * t51;
t70 = Ifges(6,6) * t52;
t69 = t17 * t49;
t65 = t48 * t49;
t64 = t49 * t51;
t26 = -t52 * mrSges(5,1) + t49 * mrSges(5,2);
t63 = mrSges(4,1) - t26;
t62 = t48 ^ 2 + t51 ^ 2;
t5 = t16 * t49 + t8 * t52;
t60 = t5 * t52 + t73;
t58 = mrSges(6,1) * t48 + mrSges(6,2) * t51;
t24 = -t52 * pkin(4) - t49 * pkin(9) - pkin(3);
t11 = t51 * t24 - t48 * t75;
t12 = t48 * t24 + t51 * t75;
t57 = -t11 * t48 + t12 * t51;
t19 = t49 * t46 + t52 * t68;
t56 = t19 * t52 + t69;
t55 = pkin(8) ^ 2;
t41 = t52 ^ 2;
t39 = t49 ^ 2;
t36 = t43 ^ 2;
t35 = t39 * t55;
t34 = Ifges(6,5) * t48;
t33 = Ifges(6,6) * t51;
t31 = t36 * t53 ^ 2;
t30 = Ifges(6,5) * t64;
t28 = Ifges(6,1) * t48 + t71;
t27 = Ifges(6,2) * t51 + t72;
t23 = -t52 * mrSges(6,1) - mrSges(6,3) * t64;
t22 = t52 * mrSges(6,2) - mrSges(6,3) * t65;
t20 = t58 * t49;
t15 = -Ifges(6,5) * t52 + (Ifges(6,1) * t51 - t72) * t49;
t14 = -t70 + (-Ifges(6,2) * t48 + t71) * t49;
t10 = t51 * t19 - t48 * t67;
t9 = -t48 * t19 - t51 * t67;
t2 = t6 * t48 + t5 * t51;
t1 = -t5 * t48 + t6 * t51;
t4 = [m(2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t79) + m(5) * (t5 ^ 2 + t78 + t79) + m(4) * (t16 ^ 2 + t8 ^ 2 + t78) + m(3) * (t47 ^ 2 + (t42 ^ 2 + t45 ^ 2) * t44 ^ 2); m(3) * t47 + m(6) * (t9 * t1 + t10 * t2 + t74) + m(5) * (t19 * t5 - t6 * t67 + t74) + m(4) * (t46 * t16 + (t50 * t8 - t53 * t6) * t43); m(3) + m(4) * (t36 * t50 ^ 2 + t46 ^ 2 + t31) + m(5) * (t19 ^ 2 + t31 + t77) + m(6) * (t10 ^ 2 + t9 ^ 2 + t77); -t8 * mrSges(4,2) + t1 * t23 + t2 * t22 + t3 * t20 - t63 * t6 + t60 * mrSges(5,3) + m(6) * (pkin(8) * t73 + t11 * t1 + t12 * t2) + m(5) * (-pkin(3) * t6 + t60 * pkin(8)); t10 * t22 + t17 * t20 + t9 * t23 + t56 * mrSges(5,3) + (-t50 * mrSges(4,2) + t63 * t53) * t43 + m(5) * (pkin(3) * t67 + t56 * pkin(8)) + m(6) * (pkin(8) * t69 + t12 * t10 + t11 * t9); -0.2e1 * pkin(3) * t26 + 0.2e1 * t11 * t23 + 0.2e1 * t12 * t22 + Ifges(4,3) + (t39 + t41) * mrSges(5,3) * t82 + m(6) * (t11 ^ 2 + t12 ^ 2 + t35) + m(5) * (pkin(3) ^ 2 + t41 * t55 + t35) + (-t30 + (Ifges(6,3) + Ifges(5,2)) * t52) * t52 + (Ifges(5,1) * t49 + 0.2e1 * Ifges(5,4) * t52 + t20 * t82 + t51 * t15 + (-t14 + t70) * t48) * t49; -t5 * mrSges(5,2) + t81 * (-t1 * t48 + t2 * t51) + t80 * t3; -t19 * mrSges(5,2) + t81 * (t10 * t51 - t9 * t48) + t80 * t17; -pkin(4) * t20 + t48 * t15 / 0.2e1 + t14 * t76 + (-pkin(8) * mrSges(5,2) - t34 / 0.2e1 - t33 / 0.2e1 + Ifges(5,6)) * t52 + t57 * mrSges(6,3) + (m(6) * t57 + t51 * t22 - t48 * t23) * pkin(9) + (-t48 * t27 / 0.2e1 + t28 * t76 + Ifges(5,5) + t80 * pkin(8)) * t49; Ifges(5,3) + t48 * t28 + t51 * t27 - 0.2e1 * pkin(4) * t25 + m(6) * (t62 * pkin(9) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t62 * pkin(9) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2); t9 * mrSges(6,1) - t10 * mrSges(6,2); t11 * mrSges(6,1) - t12 * mrSges(6,2) - Ifges(6,6) * t65 - Ifges(6,3) * t52 + t30; -t58 * pkin(9) + t33 + t34; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
