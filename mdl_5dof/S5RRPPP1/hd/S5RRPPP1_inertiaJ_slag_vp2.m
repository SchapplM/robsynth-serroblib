% Calculate joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:22
% EndTime: 2019-12-31 19:23:24
% DurationCPUTime: 0.68s
% Computational Cost: add. (737->224), mult. (1755->312), div. (0->0), fcn. (1638->6), ass. (0->89)
t87 = sin(qJ(2));
t88 = cos(qJ(2));
t84 = sin(pkin(5));
t95 = qJ(3) * t84;
t59 = -t88 * pkin(2) - t87 * t95 - pkin(1);
t86 = cos(pkin(5));
t92 = qJ(3) * t86 + pkin(7);
t66 = t92 * t87;
t85 = cos(pkin(8));
t108 = (t59 * t84 - t66 * t86) * t85;
t107 = m(4) * pkin(2);
t106 = m(5) + m(6);
t21 = t86 * t59 + t84 * t66;
t105 = m(4) * t21;
t83 = sin(pkin(8));
t104 = t83 * t84;
t103 = t83 * t86;
t102 = t84 * t85;
t101 = t84 * t88;
t100 = t85 * t86;
t99 = t86 * t88;
t98 = pkin(3) + qJ(5);
t67 = t92 * t88;
t57 = t83 * t67;
t97 = pkin(3) * t101 + t57;
t56 = pkin(2) * t103 + t85 * t95;
t51 = t83 * t99 + t87 * t85;
t25 = t51 * mrSges(6,1) + mrSges(6,3) * t101;
t64 = mrSges(6,1) * t102 + t86 * mrSges(6,2);
t65 = mrSges(5,1) * t104 + t86 * mrSges(5,2);
t96 = t87 ^ 2 + t88 ^ 2;
t8 = -t66 * t103 + t59 * t104 + t85 * t67;
t94 = -pkin(2) * t85 - pkin(3);
t50 = t87 * t83 - t85 * t99;
t18 = -t51 * mrSges(6,2) + t50 * mrSges(6,3);
t62 = mrSges(6,1) * t104 - t86 * mrSges(6,3);
t93 = -qJ(4) * t83 - pkin(2);
t31 = -t86 * qJ(4) - t56;
t28 = t51 * mrSges(5,1) - mrSges(5,2) * t101;
t27 = -t50 * mrSges(6,1) - mrSges(6,2) * t101;
t90 = -t51 * qJ(4) + t21;
t5 = qJ(4) * t101 - t8;
t74 = mrSges(5,2) * t102;
t71 = mrSges(4,2) * t104;
t69 = t83 * t95;
t63 = -mrSges(5,1) * t102 - t86 * mrSges(5,3);
t61 = -t86 * mrSges(4,2) + mrSges(4,3) * t102;
t60 = t86 * mrSges(4,1) - mrSges(4,3) * t104;
t55 = (-t83 * mrSges(6,2) - t85 * mrSges(6,3)) * t84;
t54 = -mrSges(4,1) * t102 + t71;
t53 = -mrSges(5,3) * t104 + t74;
t52 = pkin(2) * t100 - t69;
t48 = t51 * mrSges(5,3);
t46 = t51 * mrSges(4,2);
t42 = (-pkin(3) * t85 + t93) * t84;
t41 = t94 * t86 + t69;
t40 = Ifges(5,1) * t86 + (-Ifges(5,4) * t83 - Ifges(5,5) * t85) * t84;
t39 = Ifges(6,1) * t86 + (-Ifges(6,4) * t85 + Ifges(6,5) * t83) * t84;
t38 = Ifges(5,4) * t86 + (-Ifges(5,2) * t83 - Ifges(5,6) * t85) * t84;
t37 = Ifges(6,4) * t86 + (-Ifges(6,2) * t85 + Ifges(6,6) * t83) * t84;
t36 = Ifges(5,5) * t86 + (-Ifges(5,6) * t83 - Ifges(5,3) * t85) * t84;
t35 = Ifges(6,5) * t86 + (-Ifges(6,6) * t85 + Ifges(6,3) * t83) * t84;
t34 = Ifges(4,5) * t86 + (Ifges(4,1) * t83 + Ifges(4,4) * t85) * t84;
t33 = Ifges(4,6) * t86 + (Ifges(4,4) * t83 + Ifges(4,2) * t85) * t84;
t32 = Ifges(4,3) * t86 + (Ifges(4,5) * t83 + Ifges(4,6) * t85) * t84;
t30 = -mrSges(4,1) * t101 - t51 * mrSges(4,3);
t29 = mrSges(4,2) * t101 - t50 * mrSges(4,3);
t26 = t50 * mrSges(5,1) + mrSges(5,3) * t101;
t24 = (-t98 * t85 + t93) * t84;
t23 = pkin(4) * t102 - t31;
t22 = pkin(4) * t104 + t69 + (-qJ(5) + t94) * t86;
t20 = -t50 * mrSges(5,2) - t48;
t19 = t50 * mrSges(4,1) + t46;
t17 = Ifges(4,1) * t51 - Ifges(4,4) * t50 - Ifges(4,5) * t101;
t16 = Ifges(4,4) * t51 - Ifges(4,2) * t50 - Ifges(4,6) * t101;
t15 = Ifges(4,5) * t51 - Ifges(4,6) * t50 - Ifges(4,3) * t101;
t14 = -Ifges(5,1) * t101 - Ifges(5,4) * t51 + Ifges(5,5) * t50;
t13 = -Ifges(6,1) * t101 + Ifges(6,4) * t50 + Ifges(6,5) * t51;
t12 = -Ifges(5,4) * t101 - Ifges(5,2) * t51 + Ifges(5,6) * t50;
t11 = -Ifges(6,4) * t101 + Ifges(6,2) * t50 + Ifges(6,6) * t51;
t10 = -Ifges(5,5) * t101 - Ifges(5,6) * t51 + Ifges(5,3) * t50;
t9 = -Ifges(6,5) * t101 + Ifges(6,6) * t50 + Ifges(6,3) * t51;
t7 = -t57 + t108;
t6 = t97 - t108;
t4 = t50 * pkin(3) + t90;
t3 = -t50 * pkin(4) - t5;
t2 = t98 * t50 + t90;
t1 = t66 * t100 + t51 * pkin(4) + (qJ(5) * t88 - t59 * t85) * t84 + t97;
t43 = [0.2e1 * t1 * t25 + 0.2e1 * t2 * t18 + 0.2e1 * t21 * t19 + 0.2e1 * t4 * t20 + 0.2e1 * t5 * t26 + 0.2e1 * t3 * t27 + 0.2e1 * t6 * t28 + 0.2e1 * t8 * t29 + 0.2e1 * t7 * t30 + Ifges(2,3) + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t87) * t87 + 0.2e1 * t96 * pkin(7) * mrSges(3,3) + (t17 - t12 + t9) * t51 + (t10 + t11 - t16) * t50 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t87 + Ifges(3,2) * t88 + (-t13 - t14 - t15) * t84) * t88 + m(3) * (t96 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t21 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2); m(4) * (t52 * t7 + t56 * t8) + m(5) * (t31 * t5 + t42 * t4 + t41 * t6) + m(6) * (t22 * t1 + t24 * t2 + t23 * t3) + (t13 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1) * t86 + (-t38 / 0.2e1 + t34 / 0.2e1 + t35 / 0.2e1) * t51 + (t36 / 0.2e1 + t37 / 0.2e1 - t33 / 0.2e1) * t50 + (-t87 * mrSges(3,1) - t88 * mrSges(3,2)) * pkin(7) + Ifges(3,5) * t87 + Ifges(3,6) * t88 + t52 * t30 + t4 * t53 + t21 * t54 + t2 * t55 + t56 * t29 + t7 * t60 + t8 * t61 + t1 * t62 + t5 * t63 + t3 * t64 + t6 * t65 + t24 * t18 + t22 * t25 + t23 * t27 + t31 * t26 + t41 * t28 + t42 * t20 + ((-t19 - t105) * pkin(2) + (-t39 / 0.2e1 - t40 / 0.2e1 - t32 / 0.2e1) * t88 + (-t10 / 0.2e1 - t11 / 0.2e1 + t16 / 0.2e1) * t85 + (t9 / 0.2e1 - t12 / 0.2e1 + t17 / 0.2e1) * t83) * t84; 0.2e1 * t22 * t62 + 0.2e1 * t23 * t64 + 0.2e1 * t24 * t55 + 0.2e1 * t31 * t63 + 0.2e1 * t41 * t65 + 0.2e1 * t42 * t53 + 0.2e1 * t52 * t60 + 0.2e1 * t56 * t61 + Ifges(3,3) + (t39 + t40 + t32) * t86 + m(6) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(4) * (t52 ^ 2 + t56 ^ 2) + m(5) * (t31 ^ 2 + t41 ^ 2 + t42 ^ 2) + ((t84 * t107 - 0.2e1 * t54) * pkin(2) + (t33 - t36 - t37) * t85 + (t34 + t35 - t38) * t83) * t84; t46 - t48 + (-mrSges(5,2) + mrSges(4,1)) * t50 + t105 + m(5) * t4 + m(6) * t2 + t18; t71 + t74 + m(5) * t42 + m(6) * t24 + (-t107 + (-mrSges(4,1) - mrSges(6,3)) * t85 + (-mrSges(6,2) - mrSges(5,3)) * t83) * t84; m(4) + t106; m(5) * t6 + m(6) * t1 + t25 + t28; m(5) * t41 + m(6) * t22 + t62 + t65; 0; t106; m(6) * t3 + t27; m(6) * t23 + t64; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t43(1), t43(2), t43(4), t43(7), t43(11); t43(2), t43(3), t43(5), t43(8), t43(12); t43(4), t43(5), t43(6), t43(9), t43(13); t43(7), t43(8), t43(9), t43(10), t43(14); t43(11), t43(12), t43(13), t43(14), t43(15);];
Mq = res;
