% Calculate joint inertia matrix for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:49
% EndTime: 2019-03-08 20:31:51
% DurationCPUTime: 0.76s
% Computational Cost: add. (1519->196), mult. (3158->295), div. (0->0), fcn. (3664->12), ass. (0->86)
t76 = sin(qJ(6));
t80 = cos(qJ(6));
t56 = -mrSges(7,1) * t80 + mrSges(7,2) * t76;
t124 = -mrSges(6,1) + t56;
t72 = sin(pkin(12));
t74 = cos(pkin(12));
t78 = sin(qJ(4));
t82 = cos(qJ(4));
t51 = -t72 * t78 + t74 * t82;
t52 = t72 * t82 + t74 * t78;
t77 = sin(qJ(5));
t81 = cos(qJ(5));
t39 = -t81 * t51 + t52 * t77;
t40 = t51 * t77 + t52 * t81;
t25 = t39 * mrSges(6,1) + t40 * mrSges(6,2);
t41 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t55 = -t74 * mrSges(4,1) + t72 * mrSges(4,2);
t123 = -t25 - t41 - t55;
t105 = pkin(8) + qJ(3);
t97 = t105 * t74;
t98 = t105 * t72;
t43 = -t78 * t98 + t82 * t97;
t31 = pkin(9) * t51 + t43;
t42 = -t78 * t97 - t82 * t98;
t88 = -t52 * pkin(9) + t42;
t14 = t31 * t77 - t81 * t88;
t94 = mrSges(7,1) * t76 + mrSges(7,2) * t80;
t22 = t94 * t40;
t122 = m(7) * t14 + t22;
t16 = t81 * t31 + t77 * t88;
t61 = -pkin(3) * t74 - pkin(2);
t44 = -pkin(4) * t51 + t61;
t17 = pkin(5) * t39 - pkin(10) * t40 + t44;
t3 = t16 * t80 + t17 * t76;
t117 = t3 * t80;
t2 = -t16 * t76 + t17 * t80;
t108 = t76 * mrSges(7,3);
t23 = -mrSges(7,2) * t39 - t40 * t108;
t111 = t40 * t80;
t24 = mrSges(7,1) * t39 - mrSges(7,3) * t111;
t121 = m(7) * (-t2 * t76 + t117) + t80 * t23 - t76 * t24;
t120 = t14 ^ 2;
t73 = sin(pkin(6));
t79 = sin(qJ(2));
t110 = t73 * t79;
t75 = cos(pkin(6));
t45 = -t72 * t110 + t74 * t75;
t46 = t74 * t110 + t72 * t75;
t32 = t45 * t82 - t46 * t78;
t33 = t45 * t78 + t46 * t82;
t19 = -t81 * t32 + t33 * t77;
t119 = t19 ^ 2;
t69 = t74 ^ 2;
t118 = 0.2e1 * t14;
t116 = Ifges(7,4) * t76;
t115 = Ifges(7,4) * t80;
t83 = cos(qJ(2));
t109 = t73 * t83;
t21 = t32 * t77 + t33 * t81;
t11 = -t76 * t109 + t21 * t80;
t114 = t11 * t80;
t113 = t14 * t19;
t112 = t40 * t76;
t104 = Ifges(7,5) * t111 + Ifges(7,3) * t39;
t103 = Ifges(7,5) * t76 + Ifges(7,6) * t80;
t102 = t72 ^ 2 + t69;
t101 = t76 ^ 2 + t80 ^ 2;
t100 = -m(4) - m(5) - m(6);
t57 = Ifges(7,2) * t80 + t116;
t58 = Ifges(7,1) * t76 + t115;
t99 = t80 * t57 + t76 * t58 + Ifges(6,3);
t62 = pkin(4) * t77 + pkin(10);
t96 = t101 * t62;
t10 = -t80 * t109 - t21 * t76;
t93 = -t10 * t76 + t114;
t92 = -t45 * t72 + t46 * t74;
t91 = 0.2e1 * t101 * mrSges(7,3);
t90 = (mrSges(6,1) * t81 - mrSges(6,2) * t77) * pkin(4);
t89 = -t21 * mrSges(6,2) + mrSges(7,3) * t114 - t10 * t108 + t124 * t19;
t6 = Ifges(7,6) * t39 + (-Ifges(7,2) * t76 + t115) * t40;
t7 = Ifges(7,5) * t39 + (Ifges(7,1) * t80 - t116) * t40;
t87 = -t16 * mrSges(6,2) + mrSges(7,3) * t117 - t2 * t108 - t57 * t112 / 0.2e1 + t58 * t111 / 0.2e1 + Ifges(6,5) * t40 + t76 * t7 / 0.2e1 + t80 * t6 / 0.2e1 + (t103 / 0.2e1 - Ifges(6,6)) * t39 + t124 * t14;
t68 = t73 ^ 2;
t63 = -pkin(4) * t81 - pkin(5);
t60 = t68 * t83 ^ 2;
t1 = [m(2) + m(7) * (t10 ^ 2 + t11 ^ 2 + t119) + m(6) * (t21 ^ 2 + t119 + t60) + m(5) * (t32 ^ 2 + t33 ^ 2 + t60) + m(4) * (t45 ^ 2 + t46 ^ 2 + t60) + m(3) * (t68 * t79 ^ 2 + t75 ^ 2 + t60); t10 * t24 + t11 * t23 + t19 * t22 + (t19 * t40 - t21 * t39) * mrSges(6,3) + (-t32 * t52 + t33 * t51) * mrSges(5,3) + t92 * mrSges(4,3) + (-t79 * mrSges(3,2) + (mrSges(3,1) + t123) * t83) * t73 + m(7) * (t10 * t2 + t11 * t3 + t113) + m(6) * (-t44 * t109 + t16 * t21 + t113) + m(5) * (-t61 * t109 + t32 * t42 + t33 * t43) + m(4) * (pkin(2) * t109 + t92 * qJ(3)); Ifges(4,2) * t69 - 0.2e1 * pkin(2) * t55 + t22 * t118 + 0.2e1 * t2 * t24 + 0.2e1 * t3 * t23 + 0.2e1 * t44 * t25 + 0.2e1 * t61 * t41 + Ifges(3,3) + (Ifges(4,1) * t72 + 0.2e1 * Ifges(4,4) * t74) * t72 + (-0.2e1 * mrSges(5,3) * t42 + Ifges(5,1) * t52) * t52 + 0.2e1 * t102 * qJ(3) * mrSges(4,3) + (0.2e1 * mrSges(5,3) * t43 + 0.2e1 * Ifges(5,4) * t52 + Ifges(5,2) * t51) * t51 + (-0.2e1 * mrSges(6,3) * t16 + Ifges(6,2) * t39 + t104) * t39 + (mrSges(6,3) * t118 + Ifges(6,1) * t40 - t76 * t6 + t80 * t7 + (-Ifges(7,6) * t76 - (2 * Ifges(6,4))) * t39) * t40 + m(4) * (t102 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2 + t61 ^ 2) + m(6) * (t16 ^ 2 + t44 ^ 2 + t120) + m(7) * (t2 ^ 2 + t3 ^ 2 + t120); m(7) * (t10 * t80 + t11 * t76) + t100 * t109; -m(4) * pkin(2) + t76 * t23 + t80 * t24 + m(7) * (t2 * t80 + t3 * t76) + m(6) * t44 + m(5) * t61 - t123; m(7) * t101 - t100; t32 * mrSges(5,1) - t33 * mrSges(5,2) + m(7) * (t63 * t19 + t93 * t62) + m(6) * (-t19 * t81 + t21 * t77) * pkin(4) + t89; Ifges(5,6) * t51 + Ifges(5,5) * t52 + t42 * mrSges(5,1) - t43 * mrSges(5,2) + t87 + (m(6) * (-t14 * t81 + t16 * t77) + (-t39 * t77 - t40 * t81) * mrSges(6,3)) * pkin(4) + t122 * t63 + t121 * t62; 0; 0.2e1 * t63 * t56 + Ifges(5,3) + 0.2e1 * t90 + t62 * t91 + m(7) * (t101 * t62 ^ 2 + t63 ^ 2) + m(6) * (t77 ^ 2 + t81 ^ 2) * pkin(4) ^ 2 + t99; m(7) * (-pkin(5) * t19 + t93 * pkin(10)) + t89; -t122 * pkin(5) + t121 * pkin(10) + t87; 0; m(7) * (-pkin(5) * t63 + pkin(10) * t96) + (-pkin(5) + t63) * t56 + t90 + (t101 * pkin(10) + t96) * mrSges(7,3) + t99; -0.2e1 * pkin(5) * t56 + m(7) * (t101 * pkin(10) ^ 2 + pkin(5) ^ 2) + pkin(10) * t91 + t99; mrSges(7,1) * t10 - mrSges(7,2) * t11; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t112 + t104; -t56; -t94 * t62 + t103; -t94 * pkin(10) + t103; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
