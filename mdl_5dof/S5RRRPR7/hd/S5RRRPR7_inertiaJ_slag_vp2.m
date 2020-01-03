% Calculate joint inertia matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:08
% EndTime: 2019-12-31 21:16:10
% DurationCPUTime: 0.64s
% Computational Cost: add. (1234->204), mult. (2390->296), div. (0->0), fcn. (2485->8), ass. (0->82)
t122 = -pkin(7) - pkin(6);
t94 = sin(qJ(2));
t109 = t122 * t94;
t97 = cos(qJ(2));
t77 = t122 * t97;
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t50 = -t96 * t109 - t93 * t77;
t127 = t50 ^ 2;
t90 = sin(pkin(9));
t91 = cos(pkin(9));
t92 = sin(qJ(5));
t95 = cos(qJ(5));
t66 = -t92 * t90 + t95 * t91;
t67 = t95 * t90 + t92 * t91;
t41 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t126 = 0.2e1 * t41;
t125 = 0.2e1 * t50;
t83 = -t97 * pkin(2) - pkin(1);
t124 = 0.2e1 * t83;
t121 = t96 * pkin(2);
t120 = Ifges(5,4) * t90;
t119 = Ifges(5,4) * t91;
t68 = t93 * t94 - t96 * t97;
t69 = t93 * t97 + t96 * t94;
t38 = t68 * pkin(3) - t69 * qJ(4) + t83;
t52 = t93 * t109 - t96 * t77;
t16 = t90 * t38 + t91 * t52;
t118 = t16 * t91;
t117 = t69 * t90;
t116 = t69 * t91;
t115 = t90 * mrSges(5,3);
t37 = mrSges(5,1) * t117 + mrSges(5,2) * t116;
t114 = Ifges(6,5) * t67 + Ifges(6,6) * t66;
t113 = t90 ^ 2 + t91 ^ 2;
t112 = t94 ^ 2 + t97 ^ 2;
t111 = 2 * mrSges(6,3);
t30 = t67 * t69;
t31 = t66 * t69;
t110 = Ifges(6,5) * t31 - Ifges(6,6) * t30 + Ifges(6,3) * t68;
t80 = -t91 * pkin(4) - pkin(3);
t73 = -t91 * mrSges(5,1) + t90 * mrSges(5,2);
t10 = t30 * mrSges(6,1) + t31 * mrSges(6,2);
t15 = t91 * t38 - t90 * t52;
t108 = t113 * qJ(4);
t42 = Ifges(6,4) * t67 + Ifges(6,2) * t66;
t43 = Ifges(6,1) * t67 + Ifges(6,4) * t66;
t75 = Ifges(5,2) * t91 + t120;
t76 = Ifges(5,1) * t90 + t119;
t107 = t66 * t42 + t67 * t43 + t91 * t75 + t90 * t76 + Ifges(4,3);
t106 = -t15 * t90 + t118;
t39 = -t68 * mrSges(5,2) - t69 * t115;
t40 = t68 * mrSges(5,1) - mrSges(5,3) * t116;
t105 = t91 * t39 - t90 * t40;
t104 = 0.2e1 * t113 * mrSges(5,3);
t103 = (t96 * mrSges(4,1) - t93 * mrSges(4,2)) * pkin(2);
t102 = t41 + t73;
t6 = t68 * pkin(4) - pkin(8) * t116 + t15;
t9 = -pkin(8) * t117 + t16;
t2 = t95 * t6 - t92 * t9;
t21 = Ifges(5,6) * t68 + (-Ifges(5,2) * t90 + t119) * t69;
t22 = Ifges(5,5) * t68 + (Ifges(5,1) * t91 - t120) * t69;
t26 = pkin(4) * t117 + t50;
t3 = t92 * t6 + t95 * t9;
t7 = Ifges(6,4) * t31 - Ifges(6,2) * t30 + Ifges(6,6) * t68;
t8 = Ifges(6,1) * t31 - Ifges(6,4) * t30 + Ifges(6,5) * t68;
t101 = -t52 * mrSges(4,2) - t15 * t115 + t66 * t7 / 0.2e1 + t26 * t41 + mrSges(5,3) * t118 - t30 * t42 / 0.2e1 + t31 * t43 / 0.2e1 + t90 * t22 / 0.2e1 + t91 * t21 / 0.2e1 - t75 * t117 / 0.2e1 + t76 * t116 / 0.2e1 + t67 * t8 / 0.2e1 - Ifges(4,6) * t68 + Ifges(4,5) * t69 + (t73 - mrSges(4,1)) * t50 + (Ifges(5,5) * t90 + Ifges(5,6) * t91 + t114) * t68 / 0.2e1 + (-t2 * t67 + t3 * t66) * mrSges(6,3);
t85 = t91 * pkin(8);
t82 = -pkin(3) - t121;
t79 = t93 * pkin(2) + qJ(4);
t74 = t91 * qJ(4) + t85;
t72 = (-pkin(8) - qJ(4)) * t90;
t71 = t80 - t121;
t57 = t91 * t79 + t85;
t56 = (-pkin(8) - t79) * t90;
t48 = t92 * t72 + t95 * t74;
t47 = t95 * t72 - t92 * t74;
t35 = t92 * t56 + t95 * t57;
t34 = t95 * t56 - t92 * t57;
t18 = t68 * mrSges(6,1) - t31 * mrSges(6,3);
t17 = -t68 * mrSges(6,2) - t30 * mrSges(6,3);
t1 = [t97 * (Ifges(3,4) * t94 + Ifges(3,2) * t97) - 0.2e1 * pkin(1) * (-t97 * mrSges(3,1) + t94 * mrSges(3,2)) + t94 * (Ifges(3,1) * t94 + Ifges(3,4) * t97) + 0.2e1 * t16 * t39 + 0.2e1 * t15 * t40 + t37 * t125 + 0.2e1 * t26 * t10 - t30 * t7 + t31 * t8 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + Ifges(2,3) + 0.2e1 * t112 * pkin(6) * mrSges(3,3) + (mrSges(4,2) * t124 + mrSges(4,3) * t125 + Ifges(4,1) * t69 - t90 * t21 + t91 * t22) * t69 + (mrSges(4,1) * t124 - 0.2e1 * t52 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t68 + (Ifges(5,5) * t91 - Ifges(5,6) * t90 - (2 * Ifges(4,4))) * t69 + t110) * t68 + m(3) * (t112 * pkin(6) ^ 2 + pkin(1) ^ 2) + m(4) * (t52 ^ 2 + t83 ^ 2 + t127) + m(5) * (t15 ^ 2 + t16 ^ 2 + t127) + m(6) * (t2 ^ 2 + t26 ^ 2 + t3 ^ 2); (m(4) * (-t50 * t96 + t52 * t93) + (-t93 * t68 - t96 * t69) * mrSges(4,3)) * pkin(2) + m(5) * (t106 * t79 + t82 * t50) + t105 * t79 + m(6) * (t34 * t2 + t71 * t26 + t35 * t3) + (-t94 * mrSges(3,1) - t97 * mrSges(3,2)) * pkin(6) + t101 + Ifges(3,6) * t97 + Ifges(3,5) * t94 + t82 * t37 + t71 * t10 + t34 * t18 + t35 * t17; t71 * t126 + 0.2e1 * t82 * t73 + Ifges(3,3) + 0.2e1 * t103 + (-t34 * t67 + t35 * t66) * t111 + t79 * t104 + m(6) * (t34 ^ 2 + t35 ^ 2 + t71 ^ 2) + m(5) * (t113 * t79 ^ 2 + t82 ^ 2) + m(4) * (t93 ^ 2 + t96 ^ 2) * pkin(2) ^ 2 + t107; m(5) * (-pkin(3) * t50 + t106 * qJ(4)) + t105 * qJ(4) + m(6) * (t47 * t2 + t80 * t26 + t48 * t3) + t101 + t80 * t10 + t47 * t18 + t48 * t17 - pkin(3) * t37; (t82 - pkin(3)) * t73 + (t71 + t80) * t41 + t103 + m(6) * (t47 * t34 + t48 * t35 + t80 * t71) + m(5) * (-pkin(3) * t82 + t79 * t108) + ((-t34 - t47) * t67 + (t35 + t48) * t66) * mrSges(6,3) + (t113 * t79 + t108) * mrSges(5,3) + t107; -0.2e1 * pkin(3) * t73 + t80 * t126 + (-t47 * t67 + t48 * t66) * t111 + qJ(4) * t104 + m(6) * (t47 ^ 2 + t48 ^ 2 + t80 ^ 2) + m(5) * (t113 * qJ(4) ^ 2 + pkin(3) ^ 2) + t107; m(5) * t50 + m(6) * t26 + t10 + t37; m(5) * t82 + m(6) * t71 + t102; -m(5) * pkin(3) + m(6) * t80 + t102; m(5) + m(6); t2 * mrSges(6,1) - t3 * mrSges(6,2) + t110; t34 * mrSges(6,1) - t35 * mrSges(6,2) + t114; t47 * mrSges(6,1) - t48 * mrSges(6,2) + t114; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
