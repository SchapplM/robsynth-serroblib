% Calculate joint inertia matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:04:46
% EndTime: 2018-11-23 16:04:47
% DurationCPUTime: 0.85s
% Computational Cost: add. (961->225), mult. (1718->307), div. (0->0), fcn. (1518->8), ass. (0->93)
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t101 = t72 ^ 2 + t75 ^ 2;
t97 = m(6) * t101;
t119 = -m(5) - t97;
t69 = sin(pkin(10));
t54 = pkin(1) * t69 + pkin(7);
t118 = pkin(4) + t54;
t76 = cos(qJ(3));
t71 = sin(qJ(6));
t74 = cos(qJ(6));
t84 = t71 * t72 - t74 * t75;
t27 = t84 * t76;
t39 = -t71 * t75 - t74 * t72;
t28 = t39 * t76;
t10 = -t27 * mrSges(7,1) + t28 * mrSges(7,2);
t106 = t75 * mrSges(6,1);
t90 = -t72 * mrSges(6,2) + t106;
t32 = t90 * t76;
t117 = t10 + t32;
t73 = sin(qJ(3));
t65 = t73 ^ 2;
t67 = t76 ^ 2;
t116 = t65 + t67;
t115 = 0.2e1 * t116;
t114 = m(5) * pkin(3);
t113 = t101 * mrSges(6,3) - mrSges(5,2);
t112 = m(7) * pkin(5);
t111 = -t72 / 0.2e1;
t77 = -pkin(3) - pkin(8);
t110 = pkin(3) * t76;
t109 = -pkin(9) + t77;
t108 = Ifges(6,4) * t72;
t107 = Ifges(6,4) * t75;
t105 = t75 * t76;
t104 = t76 * mrSges(6,3);
t70 = cos(pkin(10));
t55 = -t70 * pkin(1) - pkin(2);
t59 = t73 * qJ(4);
t94 = t55 - t59;
t24 = t77 * t76 + t94;
t33 = t118 * t73;
t9 = t75 * t24 + t72 * t33;
t46 = mrSges(6,1) * t72 + mrSges(6,2) * t75;
t103 = t46 + mrSges(5,3);
t102 = t116 * t54 ^ 2;
t34 = t118 * t76;
t99 = t39 ^ 2 + t84 ^ 2;
t98 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t73;
t96 = t101 * t77;
t95 = -mrSges(7,1) * t84 + t39 * mrSges(7,2);
t30 = t75 * t33;
t4 = t73 * pkin(5) + t30 + (pkin(9) * t76 - t24) * t72;
t5 = -pkin(9) * t105 + t9;
t2 = t4 * t74 - t5 * t71;
t3 = t4 * t71 + t5 * t74;
t92 = -t2 * t84 - t3 * t39;
t8 = -t24 * t72 + t30;
t91 = t72 * t9 + t75 * t8;
t89 = -Ifges(6,5) * t72 - Ifges(6,6) * t75;
t44 = t109 * t72;
t45 = t109 * t75;
t15 = -t44 * t71 + t45 * t74;
t16 = t44 * t74 + t45 * t71;
t88 = -t15 * t84 - t16 * t39;
t87 = -t27 * t84 - t28 * t39;
t86 = t39 * t71 + t74 * t84;
t42 = mrSges(6,1) * t73 + t72 * t104;
t43 = -mrSges(6,2) * t73 - t75 * t104;
t85 = -t75 * t42 - t72 * t43;
t36 = Ifges(7,6) * t39;
t37 = Ifges(7,5) * t84;
t83 = t15 * mrSges(7,1) - t16 * mrSges(7,2) + t36 - t37;
t82 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t98;
t81 = (mrSges(7,1) * t74 - mrSges(7,2) * t71) * pkin(5);
t78 = qJ(4) ^ 2;
t62 = Ifges(6,5) * t75;
t61 = Ifges(6,3) * t73;
t56 = pkin(5) * t72 + qJ(4);
t48 = Ifges(6,1) * t75 - t108;
t47 = -Ifges(6,2) * t72 + t107;
t31 = t94 - t110;
t26 = Ifges(6,5) * t73 + (-Ifges(6,1) * t72 - t107) * t76;
t25 = Ifges(6,6) * t73 + (-Ifges(6,2) * t75 - t108) * t76;
t20 = pkin(5) * t105 + t34;
t18 = mrSges(7,1) * t73 - mrSges(7,3) * t28;
t17 = -mrSges(7,2) * t73 + mrSges(7,3) * t27;
t14 = -Ifges(7,1) * t84 + Ifges(7,4) * t39;
t13 = -Ifges(7,4) * t84 + Ifges(7,2) * t39;
t12 = -mrSges(7,1) * t39 - mrSges(7,2) * t84;
t7 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t73;
t6 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t73;
t1 = [0.2e1 * t20 * t10 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + t27 * t6 + t28 * t7 + 0.2e1 * t34 * t32 + 0.2e1 * t8 * t42 + 0.2e1 * t9 * t43 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t55 * mrSges(4,2) - 0.2e1 * t31 * mrSges(5,3) + t61 + (Ifges(5,2) + Ifges(4,1)) * t73 + t98) * t73 + (0.2e1 * t31 * mrSges(5,2) - 0.2e1 * t55 * mrSges(4,1) - t75 * t25 - t72 * t26 + (Ifges(5,3) + Ifges(4,2)) * t76 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t89) * t73) * t76 + m(5) * (t31 ^ 2 + t102) + m(4) * (t55 ^ 2 + t102) + m(6) * (t34 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(7) * (t2 ^ 2 + t20 ^ 2 + t3 ^ 2) + (mrSges(5,1) + mrSges(4,3)) * t54 * t115 + (0.2e1 * mrSges(3,1) * t70 - 0.2e1 * mrSges(3,2) * t69 + m(3) * (t69 ^ 2 + t70 ^ 2) * pkin(1)) * pkin(1); t28 * t17 + t27 * t18 + t85 * t76 + t117 * t73 + m(7) * (t2 * t27 + t20 * t73 + t28 * t3) + m(6) * (t34 * t73 - t91 * t76); m(3) + m(7) * (t27 ^ 2 + t28 ^ 2 + t65) + m(6) * (t101 * t67 + t65) + (m(5) / 0.2e1 + m(4) / 0.2e1) * t115; t56 * t10 - t84 * t7 / 0.2e1 + t34 * t46 + qJ(4) * t32 + t39 * t6 / 0.2e1 + t20 * t12 + t27 * t13 / 0.2e1 + t28 * t14 / 0.2e1 + t16 * t17 + t15 * t18 - t92 * mrSges(7,3) + (t26 / 0.2e1 + t77 * t42 - t8 * mrSges(6,3)) * t75 + (-t25 / 0.2e1 + t77 * t43 - t9 * mrSges(6,3)) * t72 + m(6) * (qJ(4) * t34 + t91 * t77) + m(7) * (t15 * t2 + t16 * t3 + t20 * t56) + (Ifges(6,6) * t111 + t62 / 0.2e1 - t37 / 0.2e1 + t36 / 0.2e1 - Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,1)) * t73 + (-Ifges(5,5) + Ifges(4,6) + t48 * t111 - t75 * t47 / 0.2e1 + qJ(4) * mrSges(5,1)) * t76 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t76 + (-mrSges(4,1) + mrSges(5,2) - t114) * t73) * t54; -t87 * mrSges(7,3) + (mrSges(4,1) + t113) * t76 + (-mrSges(4,2) + t12 + t103) * t73 + m(7) * (t15 * t27 + t16 * t28 + t56 * t73) + m(6) * (-t76 * t96 + t59) + m(5) * (t59 + t110); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t56 * t12 + t39 * t13 - t84 * t14 - t72 * t47 + t75 * t48 + Ifges(5,1) + Ifges(4,3) + m(7) * (t15 ^ 2 + t16 ^ 2 + t56 ^ 2) + m(6) * (t101 * t77 ^ 2 + t78) + m(5) * (pkin(3) ^ 2 + t78) - 0.2e1 * mrSges(6,3) * t96 - 0.2e1 * t88 * mrSges(7,3) + 0.2e1 * t103 * qJ(4); -t39 * t17 - t84 * t18 + (m(5) * t54 + mrSges(5,1)) * t73 + m(7) * t92 + m(6) * t91 - t85; m(7) * t87 + t119 * t76; m(7) * t88 - t99 * mrSges(7,3) + t77 * t97 - t113 - t114; m(7) * t99 - t119; t8 * mrSges(6,1) - t9 * mrSges(6,2) + t61 + t89 * t76 + (m(7) * (t2 * t74 + t3 * t71) + t71 * t17 + t74 * t18) * pkin(5) + t82; (t27 * t74 + t28 * t71) * t112 - t117; t77 * t106 + t62 + (-mrSges(6,2) * t77 - Ifges(6,6)) * t72 + (m(7) * (t15 * t74 + t16 * t71) + t86 * mrSges(7,3)) * pkin(5) + t83; -t86 * t112 + t90 + t95; Ifges(6,3) + Ifges(7,3) + m(7) * (t71 ^ 2 + t74 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t81; t82; -t10; t83; t95; Ifges(7,3) + t81; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
