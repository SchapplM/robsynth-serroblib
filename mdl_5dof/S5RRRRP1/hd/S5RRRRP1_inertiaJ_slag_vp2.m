% Calculate joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:46
% EndTime: 2019-12-05 18:44:48
% DurationCPUTime: 0.53s
% Computational Cost: add. (816->142), mult. (1530->196), div. (0->0), fcn. (1570->6), ass. (0->58)
t79 = mrSges(5,2) + mrSges(6,2);
t53 = cos(qJ(3));
t43 = pkin(2) * t53 + pkin(3);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t50 = sin(qJ(3));
t75 = pkin(2) * t50;
t29 = t52 * t43 - t49 * t75;
t27 = pkin(4) + t29;
t25 = t27 * mrSges(6,1);
t26 = t29 * mrSges(5,1);
t78 = t25 + t26;
t77 = m(6) * pkin(4);
t76 = -pkin(7) - pkin(6);
t74 = pkin(3) * t49;
t73 = pkin(3) * t52;
t51 = sin(qJ(2));
t38 = t76 * t51;
t54 = cos(qJ(2));
t39 = t76 * t54;
t21 = t53 * t38 + t39 * t50;
t34 = t50 * t54 + t51 * t53;
t11 = -pkin(8) * t34 + t21;
t22 = t50 * t38 - t53 * t39;
t33 = -t50 * t51 + t53 * t54;
t12 = pkin(8) * t33 + t22;
t6 = t49 * t11 + t52 * t12;
t18 = t33 * t52 - t34 * t49;
t30 = t43 * t49 + t52 * t75;
t72 = t18 * t30;
t71 = t18 * t49;
t69 = Ifges(5,3) + Ifges(6,3);
t68 = t51 ^ 2 + t54 ^ 2;
t67 = Ifges(4,3) + t69;
t44 = -pkin(2) * t54 - pkin(1);
t5 = t52 * t11 - t12 * t49;
t65 = t79 * t30;
t42 = pkin(4) + t73;
t41 = t42 * mrSges(6,1);
t45 = mrSges(5,1) * t73;
t64 = t41 + t45 + t69;
t63 = t79 * t74;
t19 = t33 * t49 + t34 * t52;
t2 = -qJ(5) * t19 + t5;
t62 = m(6) * t2 - t19 * mrSges(6,3);
t24 = -pkin(3) * t33 + t44;
t61 = (mrSges(4,1) * t53 - mrSges(4,2) * t50) * pkin(2);
t3 = qJ(5) * t18 + t6;
t60 = t5 * mrSges(5,1) + t2 * mrSges(6,1) - t6 * mrSges(5,2) - t3 * mrSges(6,2) + (Ifges(5,5) + Ifges(6,5)) * t19 + (Ifges(5,6) + Ifges(6,6)) * t18;
t59 = t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,5) * t34 + Ifges(4,6) * t33 + t60;
t57 = pkin(3) ^ 2;
t55 = pkin(4) * mrSges(6,1);
t46 = t49 ^ 2 * t57;
t28 = t30 ^ 2;
t23 = t30 * t74;
t13 = t19 * mrSges(6,2);
t7 = -pkin(4) * t18 + t24;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t54 + mrSges(3,2) * t51) + t51 * (Ifges(3,1) * t51 + Ifges(3,4) * t54) + t54 * (Ifges(3,4) * t51 + Ifges(3,2) * t54) + t34 * (Ifges(4,1) * t34 + Ifges(4,4) * t33) + t33 * (Ifges(4,4) * t34 + Ifges(4,2) * t33) + 0.2e1 * t44 * (-mrSges(4,1) * t33 + mrSges(4,2) * t34) + 0.2e1 * t7 * t13 + Ifges(2,3) + (0.2e1 * t24 * mrSges(5,2) - 0.2e1 * t5 * mrSges(5,3) - 0.2e1 * t2 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t19) * t19 + m(6) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + m(5) * (t24 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t21 ^ 2 + t22 ^ 2 + t44 ^ 2) + m(3) * (t68 * pkin(6) ^ 2 + pkin(1) ^ 2) + 0.2e1 * (-t21 * t34 + t22 * t33) * mrSges(4,3) + 0.2e1 * t68 * pkin(6) * mrSges(3,3) + (-0.2e1 * t24 * mrSges(5,1) - 0.2e1 * t7 * mrSges(6,1) + 0.2e1 * t6 * mrSges(5,3) + 0.2e1 * t3 * mrSges(6,3) + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t19 + (Ifges(6,2) + Ifges(5,2)) * t18) * t18; t59 + (m(4) * (t21 * t53 + t22 * t50) + (t50 * t33 - t53 * t34) * mrSges(4,3)) * pkin(2) + m(6) * (t2 * t27 + t3 * t30) + m(5) * (t29 * t5 + t30 * t6) + (-mrSges(3,1) * t51 - mrSges(3,2) * t54) * pkin(6) + (-t19 * t27 + t72) * mrSges(6,3) + (-t19 * t29 + t72) * mrSges(5,3) + Ifges(3,5) * t51 + Ifges(3,6) * t54; Ifges(3,3) + 0.2e1 * t25 + 0.2e1 * t26 + m(6) * (t27 ^ 2 + t28) + m(5) * (t29 ^ 2 + t28) + m(4) * (t50 ^ 2 + t53 ^ 2) * pkin(2) ^ 2 + t67 - 0.2e1 * t65 + 0.2e1 * t61; t59 + t62 * t42 + (mrSges(6,3) * t71 + (-t19 * t52 + t71) * mrSges(5,3) + m(6) * t3 * t49 + m(5) * (t49 * t6 + t5 * t52)) * pkin(3); Ifges(4,3) + t61 + m(6) * (t27 * t42 + t23) + m(5) * (t29 * t73 + t23) + t64 + t79 * (-t30 - t74) + t78; 0.2e1 * t41 + 0.2e1 * t45 - 0.2e1 * t63 + m(5) * (t52 ^ 2 * t57 + t46) + m(6) * (t42 ^ 2 + t46) + t67; t62 * pkin(4) + t60; t27 * t77 + t55 - t65 + t69 + t78; t42 * t77 + t55 - t63 + t64; m(6) * pkin(4) ^ 2 + 0.2e1 * t55 + t69; m(6) * t7 - t18 * mrSges(6,1) + t13; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
