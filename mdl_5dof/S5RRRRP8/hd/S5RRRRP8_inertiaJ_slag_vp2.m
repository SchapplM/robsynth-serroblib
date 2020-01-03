% Calculate joint inertia matrix for
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:40
% EndTime: 2019-12-31 21:59:43
% DurationCPUTime: 0.75s
% Computational Cost: add. (905->223), mult. (1833->309), div. (0->0), fcn. (1723->6), ass. (0->91)
t114 = 2 * pkin(6);
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t113 = (t86 * mrSges(5,1) + (-mrSges(5,2) - mrSges(6,2)) * t83) * pkin(3);
t112 = 2 * mrSges(6,1);
t111 = m(6) * pkin(4);
t110 = -pkin(8) - pkin(7);
t109 = m(6) * t83;
t88 = cos(qJ(2));
t108 = pkin(6) * t88;
t85 = sin(qJ(2));
t78 = t85 * pkin(6);
t84 = sin(qJ(3));
t107 = Ifges(4,4) * t84;
t87 = cos(qJ(3));
t106 = Ifges(4,4) * t87;
t55 = -t83 * t84 + t86 * t87;
t105 = t55 * t83;
t104 = t84 * t85;
t103 = t85 * t87;
t101 = Ifges(5,3) + Ifges(6,3);
t63 = -t88 * pkin(2) - t85 * pkin(7) - pkin(1);
t54 = t87 * t63;
t18 = -pkin(8) * t103 + t54 + (-pkin(6) * t84 - pkin(3)) * t88;
t36 = t87 * t108 + t84 * t63;
t27 = -pkin(8) * t104 + t36;
t6 = t83 * t18 + t86 * t27;
t67 = t110 * t84;
t68 = t110 * t87;
t30 = t83 * t67 - t86 * t68;
t62 = pkin(3) * t104 + t78;
t100 = t84 ^ 2 + t87 ^ 2;
t99 = -Ifges(4,3) - t101;
t46 = t55 * t85;
t5 = t86 * t18 - t83 * t27;
t2 = -t88 * pkin(4) - t46 * qJ(5) + t5;
t33 = -t88 * mrSges(6,1) - t46 * mrSges(6,3);
t98 = m(6) * t2 + t33;
t73 = -t87 * pkin(3) - pkin(2);
t56 = t83 * t87 + t86 * t84;
t45 = t56 * t85;
t14 = t45 * mrSges(6,1) + t46 * mrSges(6,2);
t21 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t29 = t86 * t67 + t83 * t68;
t96 = (-Ifges(5,5) - Ifges(6,5)) * t46 + (Ifges(5,6) + Ifges(6,6)) * t45;
t12 = -t56 * qJ(5) + t29;
t95 = m(6) * t12 - t56 * mrSges(6,3);
t94 = t84 * mrSges(4,1) + t87 * mrSges(4,2);
t3 = -t45 * qJ(5) + t6;
t93 = t5 * mrSges(5,1) + t2 * mrSges(6,1) - t6 * mrSges(5,2) - t3 * mrSges(6,2) - t96;
t13 = t55 * qJ(5) + t30;
t49 = Ifges(6,6) * t55;
t50 = Ifges(5,6) * t55;
t51 = Ifges(6,5) * t56;
t52 = Ifges(5,5) * t56;
t92 = t29 * mrSges(5,1) + t12 * mrSges(6,1) - t30 * mrSges(5,2) - t13 * mrSges(6,2) + t49 + t50 + t51 + t52;
t91 = pkin(3) ^ 2;
t90 = pkin(6) ^ 2;
t82 = t88 ^ 2;
t80 = t85 ^ 2;
t77 = t80 * t90;
t76 = t83 ^ 2 * t91;
t75 = Ifges(4,5) * t84;
t74 = Ifges(4,6) * t87;
t72 = t86 * pkin(3) + pkin(4);
t69 = Ifges(4,5) * t103;
t66 = Ifges(4,1) * t84 + t106;
t65 = Ifges(4,2) * t87 + t107;
t64 = -t87 * mrSges(4,1) + t84 * mrSges(4,2);
t61 = -t88 * mrSges(4,1) - mrSges(4,3) * t103;
t60 = t88 * mrSges(4,2) - mrSges(4,3) * t104;
t47 = t94 * t85;
t44 = -Ifges(4,5) * t88 + (Ifges(4,1) * t87 - t107) * t85;
t43 = -Ifges(4,6) * t88 + (-Ifges(4,2) * t84 + t106) * t85;
t37 = -t55 * pkin(4) + t73;
t35 = -t84 * t108 + t54;
t34 = -t88 * mrSges(5,1) - t46 * mrSges(5,3);
t32 = t88 * mrSges(5,2) - t45 * mrSges(5,3);
t31 = t88 * mrSges(6,2) - t45 * mrSges(6,3);
t26 = Ifges(5,1) * t56 + Ifges(5,4) * t55;
t25 = Ifges(6,1) * t56 + Ifges(6,4) * t55;
t24 = Ifges(5,4) * t56 + Ifges(5,2) * t55;
t23 = Ifges(6,4) * t56 + Ifges(6,2) * t55;
t22 = -t55 * mrSges(5,1) + t56 * mrSges(5,2);
t19 = t45 * pkin(4) + t62;
t15 = t45 * mrSges(5,1) + t46 * mrSges(5,2);
t10 = Ifges(5,1) * t46 - Ifges(5,4) * t45 - Ifges(5,5) * t88;
t9 = Ifges(6,1) * t46 - Ifges(6,4) * t45 - Ifges(6,5) * t88;
t8 = Ifges(5,4) * t46 - Ifges(5,2) * t45 - Ifges(5,6) * t88;
t7 = Ifges(6,4) * t46 - Ifges(6,2) * t45 - Ifges(6,6) * t88;
t1 = [0.2e1 * t19 * t14 + 0.2e1 * t62 * t15 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t31 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t34 + 0.2e1 * t35 * t61 + 0.2e1 * t36 * t60 + Ifges(2,3) + (t9 + t10) * t46 - (t7 + t8) * t45 + (t80 + t82) * mrSges(3,3) * t114 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t85 + t114 * t47 - t84 * t43 + t87 * t44) * t85 + m(3) * (pkin(1) ^ 2 + t82 * t90 + t77) + m(4) * (t35 ^ 2 + t36 ^ 2 + t77) + m(5) * (t5 ^ 2 + t6 ^ 2 + t62 ^ 2) + m(6) * (t19 ^ 2 + t2 ^ 2 + t3 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t69 + (Ifges(3,2) - t99) * t88 + (Ifges(4,6) * t84 + (2 * Ifges(3,4))) * t85 + t96) * t88; -pkin(2) * t47 + t12 * t33 + t13 * t31 + t37 * t14 + t73 * t15 + t19 * t21 + t62 * t22 + t29 * t34 + t30 * t32 + (-t75 / 0.2e1 - t74 / 0.2e1 - t51 / 0.2e1 - t49 / 0.2e1 - t52 / 0.2e1 - t50 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t88 + (t25 / 0.2e1 + t26 / 0.2e1) * t46 - (t23 / 0.2e1 + t24 / 0.2e1) * t45 + (t43 / 0.2e1 + pkin(7) * t60 + t36 * mrSges(4,3)) * t87 + (t44 / 0.2e1 - pkin(7) * t61 - t35 * mrSges(4,3)) * t84 + (Ifges(3,5) - t84 * t65 / 0.2e1 + t87 * t66 / 0.2e1 + (t64 - mrSges(3,1)) * pkin(6)) * t85 + m(4) * (-pkin(2) * t78 + (-t35 * t84 + t36 * t87) * pkin(7)) + m(5) * (t29 * t5 + t30 * t6 + t73 * t62) + m(6) * (t12 * t2 + t13 * t3 + t19 * t37) + (-t2 * mrSges(6,3) - t5 * mrSges(5,3) + t9 / 0.2e1 + t10 / 0.2e1) * t56 + (t3 * mrSges(6,3) + t6 * mrSges(5,3) + t7 / 0.2e1 + t8 / 0.2e1) * t55; -0.2e1 * pkin(2) * t64 + 0.2e1 * t37 * t21 + 0.2e1 * t73 * t22 + t87 * t65 + t84 * t66 + Ifges(3,3) + 0.2e1 * t100 * pkin(7) * mrSges(4,3) + m(6) * (t12 ^ 2 + t13 ^ 2 + t37 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t73 ^ 2) + m(4) * (t100 * pkin(7) ^ 2 + pkin(2) ^ 2) + (-0.2e1 * mrSges(5,3) * t29 - 0.2e1 * mrSges(6,3) * t12 + t25 + t26) * t56 + (0.2e1 * mrSges(5,3) * t30 + 0.2e1 * mrSges(6,3) * t13 + t23 + t24) * t55; -Ifges(4,6) * t104 + t35 * mrSges(4,1) - t36 * mrSges(4,2) + t69 + t98 * t72 + (t86 * t34 + (t31 + t32) * t83 + t3 * t109 + m(5) * (t5 * t86 + t6 * t83)) * pkin(3) + t99 * t88 + t93; t74 + t75 + t95 * t72 - t94 * pkin(7) + (mrSges(6,3) * t105 + (-t56 * t86 + t105) * mrSges(5,3) + t13 * t109 + m(5) * (t29 * t86 + t30 * t83)) * pkin(3) + t92; t72 * t112 + m(5) * (t86 ^ 2 * t91 + t76) + m(6) * (t72 ^ 2 + t76) + 0.2e1 * t113 - t99; t98 * pkin(4) - t101 * t88 + t93; t95 * pkin(4) + t92; t72 * t111 + (pkin(4) + t72) * mrSges(6,1) + t113 + t101; (t112 + t111) * pkin(4) + t101; m(6) * t19 + t14; m(6) * t37 + t21; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
