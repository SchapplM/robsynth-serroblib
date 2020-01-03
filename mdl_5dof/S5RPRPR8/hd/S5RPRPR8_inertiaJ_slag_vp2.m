% Calculate joint inertia matrix for
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:19
% EndTime: 2019-12-31 18:21:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (545->160), mult. (1100->233), div. (0->0), fcn. (987->8), ass. (0->72)
t55 = sin(pkin(8));
t46 = t55 * pkin(1) + pkin(6);
t84 = 0.2e1 * t46;
t54 = sin(pkin(9));
t59 = sin(qJ(3));
t26 = (pkin(4) * t54 + t46) * t59;
t56 = cos(pkin(9));
t73 = t56 * t59;
t74 = t54 * t59;
t27 = mrSges(5,1) * t74 + mrSges(5,2) * t73;
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t33 = t60 * t54 + t58 * t56;
t22 = t33 * t59;
t32 = -t58 * t54 + t60 * t56;
t23 = t32 * t59;
t6 = -t22 * mrSges(6,1) - t23 * mrSges(6,2);
t83 = -m(6) * t26 - t27 + t6;
t82 = m(5) * pkin(3);
t81 = -t54 / 0.2e1;
t80 = t56 / 0.2e1;
t79 = m(5) + m(6);
t61 = cos(qJ(3));
t78 = t61 * pkin(3);
t77 = Ifges(5,4) * t54;
t76 = Ifges(5,4) * t56;
t75 = t46 * t61;
t72 = t59 * mrSges(5,3);
t38 = -t56 * mrSges(5,1) + t54 * mrSges(5,2);
t71 = mrSges(4,1) - t38;
t70 = pkin(7) + qJ(4);
t69 = -Ifges(6,5) * t23 + Ifges(6,6) * t22;
t57 = cos(pkin(8));
t47 = -t57 * pkin(1) - pkin(2);
t31 = -t59 * qJ(4) + t47 - t78;
t9 = t54 * t31 + t56 * t75;
t68 = t54 ^ 2 + t56 ^ 2;
t52 = t59 ^ 2;
t53 = t61 ^ 2;
t67 = t52 + t53;
t10 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t66 = qJ(4) * t68;
t25 = t56 * t31;
t8 = -t54 * t75 + t25;
t65 = -t54 * t8 + t56 * t9;
t34 = t61 * mrSges(5,2) - t54 * t72;
t35 = -t61 * mrSges(5,1) - t56 * t72;
t64 = t56 * t34 - t54 * t35;
t48 = -t56 * pkin(4) - pkin(3);
t45 = t46 ^ 2;
t42 = t52 * t45;
t41 = Ifges(5,1) * t54 + t76;
t40 = Ifges(5,2) * t56 + t77;
t39 = t70 * t56;
t37 = t70 * t54;
t30 = Ifges(6,5) * t33;
t29 = Ifges(6,6) * t32;
t21 = -Ifges(5,5) * t61 + (Ifges(5,1) * t56 - t77) * t59;
t20 = -Ifges(5,6) * t61 + (-Ifges(5,2) * t54 + t76) * t59;
t16 = -t61 * mrSges(6,1) - t23 * mrSges(6,3);
t15 = t61 * mrSges(6,2) - t22 * mrSges(6,3);
t14 = -t58 * t37 + t60 * t39;
t13 = -t60 * t37 - t58 * t39;
t12 = Ifges(6,1) * t33 + Ifges(6,4) * t32;
t11 = Ifges(6,4) * t33 + Ifges(6,2) * t32;
t7 = -pkin(7) * t74 + t9;
t5 = -pkin(7) * t73 + t25 + (-t46 * t54 - pkin(4)) * t61;
t4 = Ifges(6,1) * t23 - Ifges(6,4) * t22 - Ifges(6,5) * t61;
t3 = Ifges(6,4) * t23 - Ifges(6,2) * t22 - Ifges(6,6) * t61;
t2 = t58 * t5 + t60 * t7;
t1 = t60 * t5 - t58 * t7;
t17 = [0.2e1 * t1 * t16 + 0.2e1 * t2 * t15 - t22 * t3 + t23 * t4 - 0.2e1 * t26 * t6 + 0.2e1 * t9 * t34 + 0.2e1 * t8 * t35 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * t47 * mrSges(4,1) + (Ifges(5,3) + Ifges(4,2) + Ifges(6,3)) * t61 + t69) * t61 + (0.2e1 * t47 * mrSges(4,2) + Ifges(4,1) * t59 - t54 * t20 + t56 * t21 + t27 * t84 + (-Ifges(5,5) * t56 + Ifges(5,6) * t54 + (2 * Ifges(4,4))) * t61) * t59 + m(6) * (t1 ^ 2 + t2 ^ 2 + t26 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2 + t42) + m(4) * (t53 * t45 + t47 ^ 2 + t42) + m(3) * (t55 ^ 2 + t57 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t57 * mrSges(3,1) - t55 * mrSges(3,2)) * pkin(1) + t67 * mrSges(4,3) * t84; -t22 * t16 + m(6) * (-t1 * t22 + t2 * t23) + t23 * t15 + (m(5) * (t65 - t75) + t64) * t59 + t83 * t61; m(3) + m(6) * (t22 ^ 2 + t23 ^ 2 + t53) + m(4) * t67 + m(5) * (t52 * t68 + t53); t54 * t21 / 0.2e1 + t20 * t80 - t48 * t6 - t22 * t11 / 0.2e1 + t23 * t12 / 0.2e1 + t26 * t10 - pkin(3) * t27 + t32 * t3 / 0.2e1 + t33 * t4 / 0.2e1 + t14 * t15 + t13 * t16 + m(6) * (t13 * t1 + t14 * t2 + t48 * t26) + (Ifges(5,5) * t81 - Ifges(5,6) * t56 / 0.2e1 - t30 / 0.2e1 - t29 / 0.2e1 + Ifges(4,6) - t46 * mrSges(4,2)) * t61 + (-t1 * t33 + t2 * t32) * mrSges(6,3) + t65 * mrSges(5,3) + (m(5) * t65 + t64) * qJ(4) + (Ifges(4,5) + t40 * t81 + t41 * t80 + (-t71 - t82) * t46) * t59; (t22 * t33 + t23 * t32) * mrSges(6,3) + (-t10 + t71) * t61 + (mrSges(5,3) * t68 - mrSges(4,2)) * t59 + m(6) * (-t13 * t22 + t14 * t23 - t48 * t61) + m(5) * (t59 * t66 + t78); -0.2e1 * pkin(3) * t38 + 0.2e1 * t48 * t10 + t32 * t11 + t33 * t12 + t56 * t40 + t54 * t41 + Ifges(4,3) + m(6) * (t13 ^ 2 + t14 ^ 2 + t48 ^ 2) + m(5) * (qJ(4) ^ 2 * t68 + pkin(3) ^ 2) + 0.2e1 * (-t13 * t33 + t14 * t32) * mrSges(6,3) + 0.2e1 * mrSges(5,3) * t66; m(5) * t59 * t46 - t83; -t79 * t61; m(6) * t48 + t10 + t38 - t82; t79; t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,3) * t61 - t69; t6; t13 * mrSges(6,1) - t14 * mrSges(6,2) + t29 + t30; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
