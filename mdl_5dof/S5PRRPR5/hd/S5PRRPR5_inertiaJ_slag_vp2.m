% Calculate joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:22
% EndTime: 2019-12-05 16:25:23
% DurationCPUTime: 0.45s
% Computational Cost: add. (528->148), mult. (1160->229), div. (0->0), fcn. (1190->10), ass. (0->64)
t47 = cos(pkin(10));
t79 = t47 * pkin(3);
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t27 = -t52 * mrSges(6,1) + t49 * mrSges(6,2);
t78 = -mrSges(5,1) + t27;
t48 = cos(pkin(5));
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t46 = sin(pkin(5));
t51 = sin(qJ(2));
t68 = t46 * t51;
t20 = t48 * t53 - t50 * t68;
t21 = t48 * t50 + t53 * t68;
t45 = sin(pkin(10));
t7 = -t47 * t20 + t45 * t21;
t77 = t7 ^ 2;
t66 = -qJ(4) - pkin(7);
t29 = t66 * t53;
t61 = t66 * t50;
t15 = -t45 * t29 - t47 * t61;
t76 = t15 ^ 2;
t44 = t53 ^ 2;
t75 = 0.2e1 * t15;
t74 = t52 / 0.2e1;
t73 = t15 * t7;
t72 = Ifges(6,4) * t49;
t71 = Ifges(6,4) * t52;
t25 = t45 * t53 + t47 * t50;
t70 = t25 * t49;
t69 = t25 * t52;
t54 = cos(qJ(2));
t67 = t46 * t54;
t24 = t45 * t50 - t47 * t53;
t65 = Ifges(6,5) * t69 + Ifges(6,3) * t24;
t64 = Ifges(6,5) * t49 + Ifges(6,6) * t52;
t63 = t49 ^ 2 + t52 ^ 2;
t62 = t50 ^ 2 + t44;
t37 = -t53 * pkin(3) - pkin(2);
t14 = t24 * mrSges(5,1) + t25 * mrSges(5,2);
t11 = t24 * pkin(4) - t25 * pkin(8) + t37;
t17 = -t47 * t29 + t45 * t61;
t1 = t52 * t11 - t49 * t17;
t2 = t49 * t11 + t52 * t17;
t60 = -t1 * t49 + t2 * t52;
t9 = t45 * t20 + t47 * t21;
t3 = -t49 * t9 - t52 * t67;
t4 = -t49 * t67 + t52 * t9;
t59 = -t3 * t49 + t4 * t52;
t58 = mrSges(6,1) * t49 + mrSges(6,2) * t52;
t57 = -t20 * t50 + t21 * t53;
t40 = t46 ^ 2;
t36 = -pkin(4) - t79;
t35 = t45 * pkin(3) + pkin(8);
t33 = t40 * t54 ^ 2;
t31 = Ifges(6,1) * t49 + t71;
t30 = Ifges(6,2) * t52 + t72;
t28 = -t53 * mrSges(4,1) + t50 * mrSges(4,2);
t13 = t24 * mrSges(6,1) - mrSges(6,3) * t69;
t12 = -t24 * mrSges(6,2) - mrSges(6,3) * t70;
t10 = t58 * t25;
t6 = Ifges(6,5) * t24 + (Ifges(6,1) * t52 - t72) * t25;
t5 = Ifges(6,6) * t24 + (-Ifges(6,2) * t49 + t71) * t25;
t8 = [m(2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t77) + m(5) * (t9 ^ 2 + t33 + t77) + m(4) * (t20 ^ 2 + t21 ^ 2 + t33) + m(3) * (t40 * t51 ^ 2 + t48 ^ 2 + t33); t7 * t10 + t4 * t12 + t3 * t13 + (-t9 * t24 + t7 * t25) * mrSges(5,3) + t57 * mrSges(4,3) + (-t51 * mrSges(3,2) + (mrSges(3,1) - t14 - t28) * t54) * t46 + m(6) * (t1 * t3 + t2 * t4 + t73) + m(5) * (t17 * t9 - t37 * t67 + t73) + m(4) * (pkin(2) * t67 + pkin(7) * t57); Ifges(4,2) * t44 - 0.2e1 * pkin(2) * t28 + 0.2e1 * t1 * t13 + t10 * t75 + 0.2e1 * t2 * t12 + 0.2e1 * t37 * t14 + Ifges(3,3) + (Ifges(4,1) * t50 + 0.2e1 * Ifges(4,4) * t53) * t50 + 0.2e1 * t62 * pkin(7) * mrSges(4,3) + (-0.2e1 * t17 * mrSges(5,3) + Ifges(5,2) * t24 + t65) * t24 + m(6) * (t1 ^ 2 + t2 ^ 2 + t76) + m(5) * (t17 ^ 2 + t37 ^ 2 + t76) + m(4) * (pkin(7) ^ 2 * t62 + pkin(2) ^ 2) + (mrSges(5,3) * t75 + Ifges(5,1) * t25 - t49 * t5 + t52 * t6 + (-Ifges(6,6) * t49 - (2 * Ifges(5,4))) * t24) * t25; t20 * mrSges(4,1) - t21 * mrSges(4,2) - t9 * mrSges(5,2) + t78 * t7 + t59 * mrSges(6,3) + m(6) * (t35 * t59 + t36 * t7) + m(5) * (t45 * t9 - t47 * t7) * pkin(3); Ifges(4,5) * t50 + Ifges(4,6) * t53 - Ifges(5,6) * t24 + t24 * t64 / 0.2e1 + t49 * t6 / 0.2e1 + t5 * t74 + t36 * t10 - t17 * mrSges(5,2) + (-t50 * mrSges(4,1) - t53 * mrSges(4,2)) * pkin(7) + t60 * mrSges(6,3) + (Ifges(5,5) + t31 * t74 - t49 * t30 / 0.2e1) * t25 + (m(5) * t17 * t45 + (-t45 * t24 - t47 * t25) * mrSges(5,3)) * pkin(3) + (-m(5) * t79 + m(6) * t36 + t78) * t15 + (m(6) * t60 + t52 * t12 - t49 * t13) * t35; 0.2e1 * t36 * t27 + t52 * t30 + t49 * t31 + Ifges(4,3) + Ifges(5,3) + m(6) * (t35 ^ 2 * t63 + t36 ^ 2) + m(5) * (t45 ^ 2 + t47 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (t47 * mrSges(5,1) - t45 * mrSges(5,2)) * pkin(3) + 0.2e1 * t63 * t35 * mrSges(6,3); m(6) * (t52 * t3 + t49 * t4) - m(5) * t67; t49 * t12 + t52 * t13 + m(6) * (t52 * t1 + t49 * t2) + m(5) * t37 + t14; 0; m(6) * t63 + m(5); t3 * mrSges(6,1) - t4 * mrSges(6,2); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,6) * t70 + t65; -t35 * t58 + t64; -t27; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
