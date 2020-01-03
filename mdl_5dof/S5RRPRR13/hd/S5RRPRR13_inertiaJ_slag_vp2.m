% Calculate joint inertia matrix for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:10
% EndTime: 2019-12-31 20:32:12
% DurationCPUTime: 0.74s
% Computational Cost: add. (1178->214), mult. (2415->314), div. (0->0), fcn. (2485->8), ass. (0->91)
t116 = 2 * pkin(6);
t89 = sin(pkin(9));
t90 = cos(pkin(9));
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t67 = t95 * t89 + t92 * t90;
t93 = sin(qJ(2));
t55 = t67 * t93;
t66 = -t92 * t89 + t95 * t90;
t56 = t66 * t93;
t115 = -Ifges(5,5) * t56 + Ifges(5,6) * t55;
t114 = -t89 / 0.2e1;
t113 = t90 / 0.2e1;
t96 = cos(qJ(2));
t112 = pkin(6) * t96;
t84 = t93 * pkin(6);
t111 = Ifges(4,4) * t89;
t110 = Ifges(4,4) * t90;
t109 = t89 * t93;
t108 = t90 * t93;
t107 = -Ifges(6,3) - Ifges(5,3);
t106 = pkin(7) + qJ(3);
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t24 = -t94 * t55 - t91 * t56;
t25 = -t91 * t55 + t94 * t56;
t105 = -Ifges(6,5) * t25 - Ifges(6,6) * t24;
t71 = -t96 * pkin(2) - t93 * qJ(3) - pkin(1);
t62 = t90 * t71;
t36 = -pkin(7) * t108 + t62 + (-pkin(6) * t89 - pkin(3)) * t96;
t48 = t90 * t112 + t89 * t71;
t42 = -pkin(7) * t109 + t48;
t16 = t92 * t36 + t95 * t42;
t72 = t106 * t89;
t74 = t106 * t90;
t44 = -t92 * t72 + t95 * t74;
t57 = mrSges(4,1) * t109 + mrSges(4,2) * t108;
t70 = pkin(3) * t109 + t84;
t104 = t89 ^ 2 + t90 ^ 2;
t81 = -t90 * pkin(3) - pkin(2);
t73 = -t90 * mrSges(4,1) + t89 * mrSges(4,2);
t28 = t55 * mrSges(5,1) + t56 * mrSges(5,2);
t39 = -t66 * mrSges(5,1) + t67 * mrSges(5,2);
t8 = -t24 * mrSges(6,1) + t25 * mrSges(6,2);
t34 = t94 * t66 - t91 * t67;
t35 = t91 * t66 + t94 * t67;
t12 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t15 = t95 * t36 - t92 * t42;
t43 = -t95 * t72 - t92 * t74;
t47 = -t89 * t112 + t62;
t103 = -t47 * t89 + t48 * t90;
t11 = -t55 * pkin(8) + t16;
t7 = -t96 * pkin(4) - t56 * pkin(8) + t15;
t2 = -t91 * t11 + t94 * t7;
t3 = t94 * t11 + t91 * t7;
t102 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t105;
t26 = -t67 * pkin(8) + t43;
t27 = t66 * pkin(8) + t44;
t10 = t91 * t26 + t94 * t27;
t32 = Ifges(6,6) * t34;
t33 = Ifges(6,5) * t35;
t9 = t94 * t26 - t91 * t27;
t101 = t9 * mrSges(6,1) - t10 * mrSges(6,2) + t32 + t33;
t100 = (t94 * mrSges(6,1) - t91 * mrSges(6,2)) * pkin(4);
t98 = pkin(6) ^ 2;
t88 = t96 ^ 2;
t87 = t93 ^ 2;
t83 = t87 * t98;
t76 = Ifges(4,1) * t89 + t110;
t75 = Ifges(4,2) * t90 + t111;
t69 = -t96 * mrSges(4,1) - mrSges(4,3) * t108;
t68 = t96 * mrSges(4,2) - mrSges(4,3) * t109;
t60 = Ifges(5,5) * t67;
t59 = Ifges(5,6) * t66;
t54 = -Ifges(4,5) * t96 + (Ifges(4,1) * t90 - t111) * t93;
t53 = -Ifges(4,6) * t96 + (-Ifges(4,2) * t89 + t110) * t93;
t49 = -t66 * pkin(4) + t81;
t46 = -t96 * mrSges(5,1) - t56 * mrSges(5,3);
t45 = t96 * mrSges(5,2) - t55 * mrSges(5,3);
t41 = Ifges(5,1) * t67 + Ifges(5,4) * t66;
t40 = Ifges(5,4) * t67 + Ifges(5,2) * t66;
t38 = t55 * pkin(4) + t70;
t23 = Ifges(5,1) * t56 - Ifges(5,4) * t55 - Ifges(5,5) * t96;
t22 = Ifges(5,4) * t56 - Ifges(5,2) * t55 - Ifges(5,6) * t96;
t18 = -t96 * mrSges(6,1) - t25 * mrSges(6,3);
t17 = t96 * mrSges(6,2) + t24 * mrSges(6,3);
t14 = Ifges(6,1) * t35 + Ifges(6,4) * t34;
t13 = Ifges(6,4) * t35 + Ifges(6,2) * t34;
t5 = Ifges(6,1) * t25 + Ifges(6,4) * t24 - Ifges(6,5) * t96;
t4 = Ifges(6,4) * t25 + Ifges(6,2) * t24 - Ifges(6,6) * t96;
t1 = [0.2e1 * t15 * t46 + 0.2e1 * t16 * t45 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t55 * t22 + t56 * t23 + t24 * t4 + t25 * t5 + 0.2e1 * t70 * t28 + 0.2e1 * t38 * t8 + 0.2e1 * t47 * t69 + 0.2e1 * t48 * t68 + Ifges(2,3) + (t87 + t88) * mrSges(3,3) * t116 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t93 + t57 * t116 - t89 * t53 + t90 * t54) * t93 + m(3) * (pkin(1) ^ 2 + t88 * t98 + t83) + m(5) * (t15 ^ 2 + t16 ^ 2 + t70 ^ 2) + m(4) * (t47 ^ 2 + t48 ^ 2 + t83) + m(6) * (t2 ^ 2 + t3 ^ 2 + t38 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2) - t107) * t96 + (-Ifges(4,5) * t90 + Ifges(4,6) * t89 + (2 * Ifges(3,4))) * t93 + t105 + t115) * t96; (t90 * t68 - t89 * t69) * qJ(3) + m(5) * (t43 * t15 + t44 * t16 + t81 * t70) + m(6) * (t10 * t3 + t9 * t2 + t49 * t38) + (-t2 * t35 + t3 * t34) * mrSges(6,3) + (-t15 * t67 + t16 * t66) * mrSges(5,3) + (t75 * t114 + t76 * t113 + Ifges(3,5) + (t73 - mrSges(3,1)) * pkin(6)) * t93 + (-pkin(6) * mrSges(3,2) + Ifges(4,5) * t114 - Ifges(4,6) * t90 / 0.2e1 - t60 / 0.2e1 - t59 / 0.2e1 - t33 / 0.2e1 - t32 / 0.2e1 + Ifges(3,6)) * t96 + m(4) * (-pkin(2) * t84 + t103 * qJ(3)) + t103 * mrSges(4,3) + t89 * t54 / 0.2e1 + t67 * t23 / 0.2e1 + t70 * t39 + t81 * t28 + t56 * t41 / 0.2e1 - pkin(2) * t57 + t66 * t22 / 0.2e1 + t44 * t45 + t43 * t46 + t49 * t8 - t55 * t40 / 0.2e1 + t34 * t4 / 0.2e1 + t35 * t5 / 0.2e1 + t38 * t12 + t10 * t17 + t9 * t18 + t24 * t13 / 0.2e1 + t25 * t14 / 0.2e1 + t53 * t113; -0.2e1 * pkin(2) * t73 + 0.2e1 * t49 * t12 + t34 * t13 + t35 * t14 + 0.2e1 * t81 * t39 + t66 * t40 + t67 * t41 + t90 * t75 + t89 * t76 + Ifges(3,3) + m(6) * (t10 ^ 2 + t49 ^ 2 + t9 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2 + t81 ^ 2) + m(4) * (t104 * qJ(3) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (t10 * t34 - t9 * t35) * mrSges(6,3) + 0.2e1 * (-t43 * t67 + t44 * t66) * mrSges(5,3) + 0.2e1 * t104 * qJ(3) * mrSges(4,3); m(4) * t84 + m(5) * t70 + m(6) * t38 + t28 + t57 + t8; -m(4) * pkin(2) + m(5) * t81 + m(6) * t49 + t12 + t39 + t73; m(4) + m(5) + m(6); t15 * mrSges(5,1) - t16 * mrSges(5,2) + t107 * t96 + (m(6) * (t2 * t94 + t3 * t91) + t91 * t17 + t94 * t18) * pkin(4) + t102 - t115; t43 * mrSges(5,1) - t44 * mrSges(5,2) + t59 + t60 + (m(6) * (t10 * t91 + t9 * t94) + (t91 * t34 - t94 * t35) * mrSges(6,3)) * pkin(4) + t101; 0; m(6) * (t91 ^ 2 + t94 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t100 - t107; -Ifges(6,3) * t96 + t102; t101; 0; Ifges(6,3) + t100; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
