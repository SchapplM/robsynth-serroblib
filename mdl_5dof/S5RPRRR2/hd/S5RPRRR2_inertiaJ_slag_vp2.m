% Calculate joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:18
% EndTime: 2019-12-05 18:11:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (964->118), mult. (1822->177), div. (0->0), fcn. (2057->8), ass. (0->54)
t46 = cos(pkin(9));
t44 = t46 ^ 2;
t45 = sin(pkin(9));
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t35 = -t49 * t45 + t52 * t46;
t36 = t52 * t45 + t49 * t46;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t23 = t51 * t35 - t48 * t36;
t24 = t48 * t35 + t51 * t36;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t12 = t50 * t23 - t47 * t24;
t73 = 0.2e1 * t12;
t72 = 0.2e1 * t23;
t71 = 0.2e1 * t35;
t70 = pkin(3) * t48;
t40 = t51 * pkin(3) + pkin(4);
t30 = t47 * t40 + t50 * t70;
t69 = t30 * mrSges(6,2);
t68 = t47 * mrSges(6,2);
t67 = Ifges(5,3) + Ifges(6,3);
t66 = pkin(6) + qJ(2);
t37 = t66 * t45;
t38 = t66 * t46;
t25 = -t52 * t37 - t49 * t38;
t18 = -t36 * pkin(7) + t25;
t26 = -t49 * t37 + t52 * t38;
t19 = t35 * pkin(7) + t26;
t8 = t48 * t18 + t51 * t19;
t65 = t45 ^ 2 + t44;
t64 = pkin(4) * t68;
t39 = -t46 * pkin(2) - pkin(1);
t13 = t47 * t23 + t50 * t24;
t63 = -t12 * mrSges(6,1) + t13 * mrSges(6,2);
t62 = -t46 * mrSges(3,1) + t45 * mrSges(3,2);
t61 = -t35 * mrSges(4,1) + t36 * mrSges(4,2);
t60 = -t23 * mrSges(5,1) + t24 * mrSges(5,2);
t7 = t51 * t18 - t48 * t19;
t29 = t50 * t40 - t47 * t70;
t28 = t29 * mrSges(6,1);
t59 = Ifges(6,3) + t28 - t69;
t4 = -t24 * pkin(8) + t7;
t5 = t23 * pkin(8) + t8;
t2 = t50 * t4 - t47 * t5;
t3 = t47 * t4 + t50 * t5;
t58 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t12;
t27 = -t35 * pkin(3) + t39;
t57 = (t51 * mrSges(5,1) - t48 * mrSges(5,2)) * pkin(3);
t56 = t7 * mrSges(5,1) - t8 * mrSges(5,2) + Ifges(5,5) * t24 + Ifges(5,6) * t23 + t58;
t41 = t50 * pkin(4) * mrSges(6,1);
t14 = -t23 * pkin(4) + t27;
t1 = [-0.2e1 * pkin(1) * t62 + Ifges(3,2) * t44 + Ifges(4,2) * t35 ^ 2 + 0.2e1 * t39 * t61 + Ifges(5,2) * t23 ^ 2 + 0.2e1 * t27 * t60 + Ifges(6,2) * t12 ^ 2 + 0.2e1 * t14 * t63 + Ifges(2,3) + t3 * mrSges(6,3) * t73 + t8 * mrSges(5,3) * t72 + t26 * mrSges(4,3) * t71 + (Ifges(3,1) * t45 + 0.2e1 * Ifges(3,4) * t46) * t45 + 0.2e1 * t65 * qJ(2) * mrSges(3,3) + (-0.2e1 * t25 * mrSges(4,3) + Ifges(4,1) * t36 + Ifges(4,4) * t71) * t36 + (-0.2e1 * t7 * mrSges(5,3) + Ifges(5,1) * t24 + Ifges(5,4) * t72) * t24 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t13 + Ifges(6,4) * t73) * t13 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t27 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t25 ^ 2 + t26 ^ 2 + t39 ^ 2) + m(3) * (t65 * qJ(2) ^ 2 + pkin(1) ^ 2); -m(3) * pkin(1) + m(4) * t39 + m(5) * t27 + m(6) * t14 + t60 + t61 + t62 + t63; m(3) + m(4) + m(5) + m(6); m(6) * (t29 * t2 + t30 * t3) + Ifges(4,5) * t36 + Ifges(4,6) * t35 - t26 * mrSges(4,2) + t25 * mrSges(4,1) + (t30 * t12 - t29 * t13) * mrSges(6,3) + (m(5) * (t48 * t8 + t51 * t7) + (t48 * t23 - t51 * t24) * mrSges(5,3)) * pkin(3) + t56; 0; -0.2e1 * t69 + Ifges(4,3) + 0.2e1 * t28 + 0.2e1 * t57 + m(6) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t48 ^ 2 + t51 ^ 2) * pkin(3) ^ 2 + t67; (m(6) * (t2 * t50 + t3 * t47) + (t47 * t12 - t50 * t13) * mrSges(6,3)) * pkin(4) + t56; 0; Ifges(5,3) + t41 + t57 + (m(6) * (t29 * t50 + t30 * t47) - t68) * pkin(4) + t59; -0.2e1 * t64 + 0.2e1 * t41 + m(6) * (t47 ^ 2 + t50 ^ 2) * pkin(4) ^ 2 + t67; t58; 0; t59; Ifges(6,3) + t41 - t64; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
