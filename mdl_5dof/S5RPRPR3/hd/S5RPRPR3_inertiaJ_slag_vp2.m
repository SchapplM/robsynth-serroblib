% Calculate joint inertia matrix for
% S5RPRPR3
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:58
% EndTime: 2020-01-03 11:36:00
% DurationCPUTime: 0.39s
% Computational Cost: add. (393->106), mult. (758->143), div. (0->0), fcn. (590->8), ass. (0->54)
t39 = sin(pkin(9));
t45 = cos(qJ(5));
t60 = t39 * t45;
t70 = Ifges(6,5) * t60;
t35 = t39 ^ 2;
t41 = cos(pkin(9));
t36 = t41 ^ 2;
t57 = t35 + t36;
t69 = t57 * mrSges(5,3);
t42 = cos(pkin(8));
t32 = pkin(1) * t42 + pkin(2);
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t40 = sin(pkin(8));
t65 = pkin(1) * t40;
t16 = t44 * t32 + t46 * t65;
t13 = qJ(4) + t16;
t68 = qJ(4) + t13;
t43 = sin(qJ(5));
t61 = t39 * t43;
t22 = mrSges(6,2) * t41 - mrSges(6,3) * t61;
t67 = 0.2e1 * t22;
t23 = -mrSges(6,1) * t41 - mrSges(6,3) * t60;
t66 = 0.2e1 * t23;
t64 = t13 * t41;
t15 = t32 * t46 - t44 * t65;
t63 = t15 * mrSges(4,1);
t62 = t16 * mrSges(4,2);
t59 = t43 * (-Ifges(6,6) * t41 + (Ifges(6,4) * t45 - Ifges(6,2) * t43) * t39);
t58 = t43 * t23;
t56 = t43 ^ 2 + t45 ^ 2;
t55 = qJ(4) * t13;
t54 = qJ(4) * t41;
t25 = -t41 * mrSges(5,1) + t39 * mrSges(5,2);
t19 = -mrSges(6,1) * t61 - mrSges(6,2) * t60;
t53 = t41 * t19 + t22 * t60;
t52 = t43 * t22 + t45 * t23 + t25;
t24 = -pkin(4) * t41 - pkin(7) * t39 - pkin(3);
t9 = -Ifges(6,6) * t61 - Ifges(6,3) * t41 + t70;
t51 = Ifges(4,3) + ((Ifges(6,1) * t45 - Ifges(6,4) * t43) * t60 + Ifges(5,1) * t39) * t39 + (0.2e1 * Ifges(5,4) * t39 + Ifges(5,2) * t41 - t70 - t9) * t41;
t50 = -0.2e1 * t39 * t19 + 0.2e1 * t69;
t49 = -t39 * t59 + t51;
t47 = qJ(4) ^ 2;
t33 = t35 * t47;
t14 = -pkin(3) - t15;
t11 = t13 ^ 2;
t8 = t35 * t11;
t7 = t24 * t43 + t45 * t54;
t6 = t24 * t45 - t43 * t54;
t5 = t35 * t55;
t3 = -t15 + t24;
t2 = t3 * t43 + t45 * t64;
t1 = t3 * t45 - t43 * t64;
t4 = [t49 + t50 * t13 + m(5) * (t11 * t36 + t14 ^ 2 + t8) + m(6) * (t1 ^ 2 + t2 ^ 2 + t8) + m(4) * (t15 ^ 2 + t16 ^ 2) + 0.2e1 * t63 - 0.2e1 * t62 + t2 * t67 + t1 * t66 + 0.2e1 * t14 * t25 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t42 * mrSges(3,1) - 0.2e1 * t40 * mrSges(3,2) + m(3) * (t40 ^ 2 + t42 ^ 2) * pkin(1)) * pkin(1); (m(6) * (-t1 * t43 + t2 * t45 - t64) - t58) * t39 + t53; m(3) + m(4) + m(5) * t57 + m(6) * (t35 * t56 + t36); t63 - t62 + (t14 - pkin(3)) * t25 + (t1 + t6) * t23 + (t2 + t7) * t22 + m(5) * (-pkin(3) * t14 + t36 * t55 + t5) + m(6) * (t1 * t6 + t2 * t7 + t5) + (-t19 * t68 - t59) * t39 + t51 + t68 * t69; (m(6) * (-t43 * t6 + t45 * t7 - t54) - t58) * t39 + t53; -0.2e1 * pkin(3) * t25 + t7 * t67 + t6 * t66 + m(6) * (t6 ^ 2 + t7 ^ 2 + t33) + m(5) * (pkin(3) ^ 2 + t36 * t47 + t33) + t50 * qJ(4) + t49; m(5) * t14 + m(6) * (t1 * t45 + t2 * t43) + t52; 0; m(6) * (t43 * t7 + t45 * t6) - m(5) * pkin(3) + t52; m(6) * t56 + m(5); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t9; t19; mrSges(6,1) * t6 - mrSges(6,2) * t7 + t9; mrSges(6,1) * t45 - mrSges(6,2) * t43; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
