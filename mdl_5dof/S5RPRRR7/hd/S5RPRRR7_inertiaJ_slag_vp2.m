% Calculate joint inertia matrix for
% S5RPRRR7
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:03:14
% DurationCPUTime: 0.60s
% Computational Cost: add. (653->176), mult. (1307->260), div. (0->0), fcn. (1170->8), ass. (0->77)
t59 = sin(qJ(3));
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t70 = t58 * mrSges(5,1) + t61 * mrSges(5,2);
t29 = t70 * t59;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t34 = t57 * t61 + t60 * t58;
t24 = t34 * t59;
t33 = -t57 * t58 + t60 * t61;
t25 = t33 * t59;
t7 = -t24 * mrSges(6,1) - t25 * mrSges(6,2);
t88 = -t29 + t7;
t55 = sin(pkin(9));
t46 = t55 * pkin(1) + pkin(6);
t87 = 0.2e1 * t46;
t86 = t61 / 0.2e1;
t85 = -pkin(8) - pkin(7);
t84 = t59 * pkin(7);
t62 = cos(qJ(3));
t83 = t62 * pkin(3);
t82 = Ifges(5,4) * t58;
t81 = Ifges(5,4) * t61;
t80 = t46 * t62;
t79 = t58 * t59;
t78 = t59 * t61;
t38 = -t61 * mrSges(5,1) + t58 * mrSges(5,2);
t77 = mrSges(4,1) - t38;
t76 = -Ifges(6,3) - Ifges(5,3);
t75 = -Ifges(6,5) * t25 + Ifges(6,6) * t24;
t56 = cos(pkin(9));
t47 = -t56 * pkin(1) - pkin(2);
t30 = t47 - t83 - t84;
t10 = t58 * t30 + t61 * t80;
t74 = t58 ^ 2 + t61 ^ 2;
t52 = t59 ^ 2;
t54 = t62 ^ 2;
t73 = t52 + t54;
t72 = t74 * mrSges(5,3);
t27 = t61 * t30;
t9 = -t58 * t80 + t27;
t71 = t10 * t61 - t58 * t9;
t35 = t62 * mrSges(5,2) - mrSges(5,3) * t79;
t36 = -t62 * mrSges(5,1) - mrSges(5,3) * t78;
t69 = t61 * t35 - t58 * t36;
t6 = -pkin(8) * t78 + t27 + (-t46 * t58 - pkin(4)) * t62;
t8 = -pkin(8) * t79 + t10;
t2 = -t57 * t8 + t60 * t6;
t3 = t57 * t6 + t60 * t8;
t68 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t75;
t42 = t85 * t58;
t43 = t85 * t61;
t15 = t60 * t42 + t57 * t43;
t16 = t57 * t42 - t60 * t43;
t31 = Ifges(6,6) * t33;
t32 = Ifges(6,5) * t34;
t67 = t15 * mrSges(6,1) - t16 * mrSges(6,2) + t31 + t32;
t66 = (t60 * mrSges(6,1) - t57 * mrSges(6,2)) * pkin(4);
t50 = Ifges(5,5) * t58;
t49 = Ifges(5,6) * t61;
t48 = -t61 * pkin(4) - pkin(3);
t45 = t46 ^ 2;
t44 = Ifges(5,5) * t78;
t41 = t52 * t45;
t40 = Ifges(5,1) * t58 + t81;
t39 = Ifges(5,2) * t61 + t82;
t28 = (pkin(4) * t58 + t46) * t59;
t23 = -Ifges(5,5) * t62 + (Ifges(5,1) * t61 - t82) * t59;
t22 = -Ifges(5,6) * t62 + (-Ifges(5,2) * t58 + t81) * t59;
t18 = -t62 * mrSges(6,1) - t25 * mrSges(6,3);
t17 = t62 * mrSges(6,2) - t24 * mrSges(6,3);
t13 = Ifges(6,1) * t34 + Ifges(6,4) * t33;
t12 = Ifges(6,4) * t34 + Ifges(6,2) * t33;
t11 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t5 = Ifges(6,1) * t25 - Ifges(6,4) * t24 - Ifges(6,5) * t62;
t4 = Ifges(6,4) * t25 - Ifges(6,2) * t24 - Ifges(6,6) * t62;
t1 = [0.2e1 * t10 * t35 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t24 * t4 + t25 * t5 - 0.2e1 * t28 * t7 + 0.2e1 * t9 * t36 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t47 * mrSges(4,2) + Ifges(4,1) * t59 - t58 * t22 + t61 * t23 + t29 * t87) * t59 + (-0.2e1 * t47 * mrSges(4,1) - t44 + (Ifges(4,2) - t76) * t62 + (Ifges(5,6) * t58 + (2 * Ifges(4,4))) * t59 + t75) * t62 + m(6) * (t2 ^ 2 + t28 ^ 2 + t3 ^ 2) + m(5) * (t10 ^ 2 + t9 ^ 2 + t41) + m(4) * (t54 * t45 + t47 ^ 2 + t41) + m(3) * (t55 ^ 2 + t56 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t56 * mrSges(3,1) - t55 * mrSges(3,2)) * pkin(1) + t73 * mrSges(4,3) * t87; m(6) * (-t2 * t24 + t3 * t25) - t24 * t18 + t25 * t17 + (m(5) * (t71 - t80) + t69) * t59 + (-m(6) * t28 + t88) * t62; m(3) + m(6) * (t24 ^ 2 + t25 ^ 2 + t54) + m(5) * (t52 * t74 + t54) + m(4) * t73; t58 * t23 / 0.2e1 + t22 * t86 - t48 * t7 + t33 * t4 / 0.2e1 + t34 * t5 / 0.2e1 - t24 * t12 / 0.2e1 + t25 * t13 / 0.2e1 + t28 * t11 - pkin(3) * t29 + t16 * t17 + t15 * t18 + m(6) * (t15 * t2 + t16 * t3 + t48 * t28) + (-t50 / 0.2e1 - t49 / 0.2e1 - t32 / 0.2e1 - t31 / 0.2e1 + Ifges(4,6) - t46 * mrSges(4,2)) * t62 + (-t2 * t34 + t3 * t33) * mrSges(6,3) + t71 * mrSges(5,3) + (m(5) * t71 + t69) * pkin(7) + (Ifges(4,5) - t58 * t39 / 0.2e1 + t40 * t86 + (-m(5) * pkin(3) - t77) * t46) * t59; (t24 * t34 + t25 * t33) * mrSges(6,3) + (-t11 + t77) * t62 + (-mrSges(4,2) + t72) * t59 + m(6) * (-t15 * t24 + t16 * t25 - t48 * t62) + m(5) * (t74 * t84 + t83); -0.2e1 * pkin(3) * t38 + 0.2e1 * t48 * t11 + t33 * t12 + t34 * t13 + t61 * t39 + t58 * t40 + Ifges(4,3) + m(6) * (t15 ^ 2 + t16 ^ 2 + t48 ^ 2) + m(5) * (pkin(7) ^ 2 * t74 + pkin(3) ^ 2) + 0.2e1 * (-t15 * t34 + t16 * t33) * mrSges(6,3) + 0.2e1 * pkin(7) * t72; -Ifges(5,6) * t79 + t9 * mrSges(5,1) - t10 * mrSges(5,2) + t44 + t76 * t62 + (m(6) * (t2 * t60 + t3 * t57) + t57 * t17 + t60 * t18) * pkin(4) + t68; m(6) * (-t24 * t60 + t25 * t57) * pkin(4) + t88; t49 + t50 - t70 * pkin(7) + (m(6) * (t15 * t60 + t16 * t57) + (t57 * t33 - t60 * t34) * mrSges(6,3)) * pkin(4) + t67; m(6) * (t57 ^ 2 + t60 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t66 - t76; -Ifges(6,3) * t62 + t68; t7; t67; Ifges(6,3) + t66; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
