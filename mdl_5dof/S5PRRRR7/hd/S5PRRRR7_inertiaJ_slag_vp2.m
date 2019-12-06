% Calculate joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:49
% EndTime: 2019-12-05 17:11:50
% DurationCPUTime: 0.33s
% Computational Cost: add. (550->130), mult. (1194->198), div. (0->0), fcn. (1197->8), ass. (0->53)
t48 = cos(qJ(3));
t40 = t48 ^ 2;
t67 = -pkin(7) - pkin(6);
t43 = sin(qJ(4));
t66 = pkin(3) * t43;
t47 = cos(qJ(4));
t35 = t47 * pkin(3) + pkin(4);
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t24 = t42 * t35 + t46 * t66;
t65 = t24 * mrSges(6,2);
t64 = t42 * mrSges(6,2);
t63 = Ifges(5,3) + Ifges(6,3);
t44 = sin(qJ(3));
t33 = t67 * t44;
t34 = t67 * t48;
t17 = t43 * t33 - t47 * t34;
t62 = t44 ^ 2 + t40;
t61 = pkin(4) * t64;
t28 = t43 * t48 + t47 * t44;
t45 = sin(qJ(2));
t21 = t28 * t45;
t27 = -t43 * t44 + t47 * t48;
t22 = t27 * t45;
t6 = -t46 * t21 - t42 * t22;
t7 = -t42 * t21 + t46 * t22;
t60 = t6 * mrSges(6,1) - t7 * mrSges(6,2);
t36 = -t48 * pkin(3) - pkin(2);
t59 = t62 * mrSges(4,3);
t16 = t47 * t33 + t43 * t34;
t23 = t46 * t35 - t42 * t66;
t20 = t23 * mrSges(6,1);
t58 = Ifges(6,3) + t20 - t65;
t57 = -t44 * mrSges(4,1) - t48 * mrSges(4,2);
t12 = t46 * t27 - t42 * t28;
t13 = t42 * t27 + t46 * t28;
t8 = -t28 * pkin(8) + t16;
t9 = t27 * pkin(8) + t17;
t2 = -t42 * t9 + t46 * t8;
t3 = t42 * t8 + t46 * t9;
t56 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t12;
t55 = (t47 * mrSges(5,1) - t43 * mrSges(5,2)) * pkin(3);
t54 = -t21 * mrSges(5,1) - t22 * mrSges(5,2) + t60;
t53 = t16 * mrSges(5,1) - t17 * mrSges(5,2) + Ifges(5,5) * t28 + Ifges(5,6) * t27 + t56;
t49 = cos(qJ(2));
t41 = t49 ^ 2;
t39 = t45 ^ 2;
t37 = t46 * pkin(4) * mrSges(6,1);
t32 = -t48 * mrSges(4,1) + t44 * mrSges(4,2);
t18 = -t27 * pkin(4) + t36;
t14 = -t27 * mrSges(5,1) + t28 * mrSges(5,2);
t4 = -t12 * mrSges(6,1) + t13 * mrSges(6,2);
t1 = [m(2) + m(6) * (t6 ^ 2 + t7 ^ 2 + t41) + m(5) * (t21 ^ 2 + t22 ^ 2 + t41) + m(4) * (t62 * t39 + t41) + m(3) * (t39 + t41); (t7 * t12 - t6 * t13) * mrSges(6,3) + (t21 * t28 + t22 * t27) * mrSges(5,3) + (-mrSges(3,2) + t59) * t45 + (mrSges(3,1) - t14 - t32 - t4) * t49 + m(6) * (-t18 * t49 + t2 * t6 + t3 * t7) + m(5) * (-t16 * t21 + t17 * t22 - t36 * t49) + m(4) * (t62 * t45 * pkin(6) + t49 * pkin(2)); Ifges(4,2) * t40 - 0.2e1 * pkin(2) * t32 + 0.2e1 * t36 * t14 + 0.2e1 * t18 * t4 + Ifges(3,3) + (Ifges(4,1) * t44 + 0.2e1 * Ifges(4,4) * t48) * t44 + (-0.2e1 * t16 * mrSges(5,3) + Ifges(5,1) * t28) * t28 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t13) * t13 + 0.2e1 * pkin(6) * t59 + (0.2e1 * t17 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t28 + Ifges(5,2) * t27) * t27 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t13 + Ifges(6,2) * t12) * t12 + m(6) * (t18 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t36 ^ 2) + m(4) * (t62 * pkin(6) ^ 2 + pkin(2) ^ 2); t57 * t45 + m(6) * (t23 * t6 + t24 * t7) + m(5) * (-t21 * t47 + t22 * t43) * pkin(3) + t54; m(6) * (t23 * t2 + t24 * t3) + Ifges(4,5) * t44 + Ifges(4,6) * t48 + t57 * pkin(6) + (t24 * t12 - t23 * t13) * mrSges(6,3) + (m(5) * (t16 * t47 + t17 * t43) + (t43 * t27 - t47 * t28) * mrSges(5,3)) * pkin(3) + t53; -0.2e1 * t65 + Ifges(4,3) + 0.2e1 * t20 + 0.2e1 * t55 + m(6) * (t23 ^ 2 + t24 ^ 2) + m(5) * (t43 ^ 2 + t47 ^ 2) * pkin(3) ^ 2 + t63; m(6) * (t42 * t7 + t46 * t6) * pkin(4) + t54; (m(6) * (t2 * t46 + t3 * t42) + (t42 * t12 - t46 * t13) * mrSges(6,3)) * pkin(4) + t53; Ifges(5,3) + t37 + t55 + (m(6) * (t23 * t46 + t24 * t42) - t64) * pkin(4) + t58; -0.2e1 * t61 + 0.2e1 * t37 + m(6) * (t42 ^ 2 + t46 ^ 2) * pkin(4) ^ 2 + t63; t60; t56; t58; Ifges(6,3) + t37 - t61; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
