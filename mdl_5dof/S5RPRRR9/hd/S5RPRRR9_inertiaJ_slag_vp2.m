% Calculate joint inertia matrix for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:19
% EndTime: 2019-12-31 19:07:21
% DurationCPUTime: 0.44s
% Computational Cost: add. (1015->139), mult. (1938->207), div. (0->0), fcn. (2144->8), ass. (0->63)
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t55 = cos(pkin(9));
t81 = pkin(6) + qJ(2);
t74 = t81 * t55;
t54 = sin(pkin(9));
t75 = t81 * t54;
t30 = -t58 * t75 + t61 * t74;
t36 = -t58 * t54 + t61 * t55;
t21 = t36 * pkin(7) + t30;
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t29 = -t58 * t74 - t61 * t75;
t37 = t61 * t54 + t58 * t55;
t66 = -t37 * pkin(7) + t29;
t10 = t57 * t21 - t60 * t66;
t28 = t57 * t36 + t60 * t37;
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t69 = mrSges(6,1) * t56 + mrSges(6,2) * t59;
t14 = t69 * t28;
t95 = m(6) * t10 + t14;
t27 = -t60 * t36 + t57 * t37;
t84 = t56 * mrSges(6,3);
t15 = -t27 * mrSges(6,2) - t28 * t84;
t85 = t28 * t59;
t16 = t27 * mrSges(6,1) - mrSges(6,3) * t85;
t12 = t60 * t21 + t57 * t66;
t44 = -t55 * pkin(2) - pkin(1);
t31 = -t36 * pkin(3) + t44;
t13 = t27 * pkin(4) - t28 * pkin(8) + t31;
t2 = -t56 * t12 + t59 * t13;
t3 = t59 * t12 + t56 * t13;
t89 = t3 * t59;
t94 = m(6) * (-t2 * t56 + t89) + t59 * t15 - t56 * t16;
t93 = t10 ^ 2;
t51 = t55 ^ 2;
t92 = 0.2e1 * t10;
t91 = 0.2e1 * t31;
t90 = 0.2e1 * t36;
t88 = Ifges(6,4) * t56;
t87 = Ifges(6,4) * t59;
t86 = t28 * t56;
t80 = Ifges(6,5) * t85 + Ifges(6,3) * t27;
t79 = Ifges(6,5) * t56 + Ifges(6,6) * t59;
t78 = t54 ^ 2 + t51;
t77 = t56 ^ 2 + t59 ^ 2;
t41 = Ifges(6,2) * t59 + t88;
t42 = Ifges(6,1) * t56 + t87;
t76 = t59 * t41 + t56 * t42 + Ifges(5,3);
t45 = t57 * pkin(3) + pkin(8);
t73 = t77 * t45;
t72 = -t55 * mrSges(3,1) + t54 * mrSges(3,2);
t71 = -t36 * mrSges(4,1) + t37 * mrSges(4,2);
t68 = 0.2e1 * mrSges(6,3) * t77;
t67 = (t60 * mrSges(5,1) - t57 * mrSges(5,2)) * pkin(3);
t40 = -t59 * mrSges(6,1) + t56 * mrSges(6,2);
t6 = Ifges(6,6) * t27 + (-Ifges(6,2) * t56 + t87) * t28;
t7 = Ifges(6,5) * t27 + (Ifges(6,1) * t59 - t88) * t28;
t65 = -t12 * mrSges(5,2) + mrSges(6,3) * t89 - t2 * t84 + t56 * t7 / 0.2e1 - t41 * t86 / 0.2e1 + t42 * t85 / 0.2e1 + Ifges(5,5) * t28 + t59 * t6 / 0.2e1 + (t79 / 0.2e1 - Ifges(5,6)) * t27 + (t40 - mrSges(5,1)) * t10;
t46 = -t60 * pkin(3) - pkin(4);
t23 = t28 * mrSges(5,2);
t1 = [-0.2e1 * pkin(1) * t72 + Ifges(3,2) * t51 + 0.2e1 * t44 * t71 + Ifges(4,2) * t36 ^ 2 + t23 * t91 + t14 * t92 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + t30 * mrSges(4,3) * t90 + (Ifges(3,1) * t54 + 0.2e1 * Ifges(3,4) * t55) * t54 + (mrSges(5,1) * t91 - 0.2e1 * t12 * mrSges(5,3) + Ifges(5,2) * t27 + t80) * t27 + 0.2e1 * t78 * qJ(2) * mrSges(3,3) + (-0.2e1 * t29 * mrSges(4,3) + Ifges(4,1) * t37 + Ifges(4,4) * t90) * t37 + (mrSges(5,3) * t92 + Ifges(5,1) * t28 - t56 * t6 + t59 * t7 + (-Ifges(6,6) * t56 - (2 * Ifges(5,4))) * t27) * t28 + m(6) * (t2 ^ 2 + t3 ^ 2 + t93) + m(5) * (t12 ^ 2 + t31 ^ 2 + t93) + m(4) * (t29 ^ 2 + t30 ^ 2 + t44 ^ 2) + m(3) * (qJ(2) ^ 2 * t78 + pkin(1) ^ 2); -m(3) * pkin(1) + t27 * mrSges(5,1) + t56 * t15 + t59 * t16 + t23 + m(6) * (t59 * t2 + t56 * t3) + m(5) * t31 + m(4) * t44 + t71 + t72; m(6) * t77 + m(3) + m(4) + m(5); t65 + (m(5) * (-t10 * t60 + t12 * t57) + (-t57 * t27 - t60 * t28) * mrSges(5,3)) * pkin(3) + Ifges(4,6) * t36 + Ifges(4,5) * t37 + t29 * mrSges(4,1) - t30 * mrSges(4,2) + t95 * t46 + t94 * t45; 0; 0.2e1 * t46 * t40 + Ifges(4,3) + 0.2e1 * t67 + t45 * t68 + m(6) * (t45 ^ 2 * t77 + t46 ^ 2) + m(5) * (t57 ^ 2 + t60 ^ 2) * pkin(3) ^ 2 + t76; -t95 * pkin(4) + t94 * pkin(8) + t65; 0; m(6) * (-pkin(4) * t46 + pkin(8) * t73) + (t46 - pkin(4)) * t40 + t67 + (pkin(8) * t77 + t73) * mrSges(6,3) + t76; -0.2e1 * pkin(4) * t40 + m(6) * (pkin(8) ^ 2 * t77 + pkin(4) ^ 2) + pkin(8) * t68 + t76; t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t86 + t80; -t40; -t45 * t69 + t79; -pkin(8) * t69 + t79; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
