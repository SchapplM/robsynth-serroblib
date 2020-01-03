% Calculate joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:58
% EndTime: 2019-12-31 19:34:59
% DurationCPUTime: 0.37s
% Computational Cost: add. (538->139), mult. (1003->189), div. (0->0), fcn. (959->6), ass. (0->57)
t41 = sin(pkin(8));
t42 = cos(pkin(8));
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t21 = t41 * t44 - t42 * t46;
t22 = t41 * t46 + t42 * t44;
t35 = -t46 * pkin(2) - pkin(1);
t49 = -t22 * qJ(4) + t35;
t9 = t21 * pkin(3) + t49;
t71 = -0.2e1 * t9;
t68 = t41 * pkin(2);
t31 = qJ(4) + t68;
t70 = t31 ^ 2;
t69 = 0.2e1 * t35;
t67 = t42 * pkin(2);
t43 = sin(qJ(5));
t45 = cos(qJ(5));
t56 = t43 ^ 2 + t45 ^ 2;
t24 = m(6) * t56;
t66 = m(5) + t24;
t65 = Ifges(6,4) * t43;
t64 = Ifges(6,4) * t45;
t63 = Ifges(6,6) * t43;
t62 = t21 * t43;
t61 = t21 * t45;
t60 = mrSges(5,1) + mrSges(4,3);
t59 = mrSges(5,2) - mrSges(4,1);
t58 = -qJ(3) - pkin(6);
t57 = t44 ^ 2 + t46 ^ 2;
t55 = Ifges(6,5) * t62 + Ifges(6,6) * t61 + Ifges(6,3) * t22;
t26 = t58 * t46;
t53 = t58 * t44;
t12 = -t41 * t26 - t42 * t53;
t14 = -t42 * t26 + t41 * t53;
t54 = t12 ^ 2 + t14 ^ 2;
t34 = -pkin(3) - t67;
t52 = t56 * mrSges(6,3);
t5 = (pkin(3) + pkin(7)) * t21 + t49;
t6 = t22 * pkin(4) + t12;
t1 = -t43 * t5 + t45 * t6;
t2 = t43 * t6 + t45 * t5;
t51 = t45 * t1 + t43 * t2;
t50 = t45 * mrSges(6,1) - t43 * mrSges(6,2);
t36 = Ifges(6,5) * t45;
t30 = -pkin(7) + t34;
t28 = Ifges(6,1) * t45 - t65;
t27 = -Ifges(6,2) * t43 + t64;
t25 = t43 * mrSges(6,1) + t45 * mrSges(6,2);
t19 = t22 * mrSges(5,3);
t18 = t22 * mrSges(4,2);
t11 = -t22 * mrSges(6,2) + mrSges(6,3) * t61;
t10 = t22 * mrSges(6,1) - mrSges(6,3) * t62;
t8 = t50 * t21;
t7 = -t21 * pkin(4) + t14;
t4 = Ifges(6,5) * t22 + (Ifges(6,1) * t43 + t64) * t21;
t3 = Ifges(6,6) * t22 + (Ifges(6,2) * t45 + t65) * t21;
t13 = [t46 * (Ifges(3,4) * t44 + Ifges(3,2) * t46) - 0.2e1 * pkin(1) * (-t46 * mrSges(3,1) + t44 * mrSges(3,2)) + t44 * (Ifges(3,1) * t44 + Ifges(3,4) * t46) + t18 * t69 + t19 * t71 - 0.2e1 * t7 * t8 + 0.2e1 * t1 * t10 + 0.2e1 * t2 * t11 + Ifges(2,3) + 0.2e1 * t57 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t69 + mrSges(5,2) * t71 + t45 * t3 + t43 * t4 + (Ifges(4,2) + Ifges(5,3)) * t21 - 0.2e1 * t60 * t14) * t21 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(5) * (t9 ^ 2 + t54) + m(4) * (t35 ^ 2 + t54) + m(3) * (t57 * pkin(6) ^ 2 + pkin(1) ^ 2) + ((Ifges(4,1) + Ifges(5,2)) * t22 + t55 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t21 + 0.2e1 * t60 * t12) * t22; Ifges(3,5) * t44 + Ifges(3,6) * t46 + t7 * t25 - t31 * t8 + (-mrSges(4,2) + mrSges(5,3)) * t14 + t59 * t12 + (-t44 * mrSges(3,1) - t46 * mrSges(3,2)) * pkin(6) + (t4 / 0.2e1 - t1 * mrSges(6,3) + t30 * t10) * t45 + (-t3 / 0.2e1 - t2 * mrSges(6,3) + t30 * t11) * t43 + m(5) * (t34 * t12 + t31 * t14) + m(6) * (t51 * t30 + t31 * t7) + m(4) * (-t12 * t42 + t14 * t41) * pkin(2) + (-t63 / 0.2e1 + t36 / 0.2e1 + Ifges(4,5) - Ifges(5,4) - mrSges(4,3) * t67 + t34 * mrSges(5,1)) * t22 + (Ifges(5,5) - Ifges(4,6) - mrSges(4,3) * t68 - t31 * mrSges(5,1) + t43 * t28 / 0.2e1 + t45 * t27 / 0.2e1) * t21; 0.2e1 * t34 * mrSges(5,2) - t43 * t27 + t45 * t28 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + m(5) * (t34 ^ 2 + t70) + m(6) * (t56 * t30 ^ 2 + t70) + m(4) * (t41 ^ 2 + t42 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t25 + mrSges(5,3)) * t31 + 0.2e1 * (t42 * mrSges(4,1) - t41 * mrSges(4,2)) * pkin(2) - 0.2e1 * t30 * t52; -t43 * t10 + t45 * t11 + t18 - t19 - t59 * t21 + m(6) * (-t43 * t1 + t45 * t2) + m(5) * t9 + m(4) * t35; 0; m(4) + t66; m(5) * t12 + m(6) * t51 + t22 * mrSges(5,1) + t45 * t10 + t43 * t11; m(5) * t34 + t30 * t24 + mrSges(5,2) - t52; 0; t66; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t55; t50 * t30 + t36 - t63; -t25; t50; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
