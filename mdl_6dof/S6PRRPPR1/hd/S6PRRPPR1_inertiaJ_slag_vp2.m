% Calculate joint inertia matrix for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:35
% EndTime: 2019-03-08 20:58:37
% DurationCPUTime: 0.80s
% Computational Cost: add. (1362->237), mult. (2865->357), div. (0->0), fcn. (3187->12), ass. (0->90)
t79 = sin(pkin(11));
t82 = cos(pkin(11));
t85 = sin(qJ(3));
t88 = cos(qJ(3));
t57 = t79 * t88 + t82 * t85;
t81 = cos(pkin(12));
t104 = t57 * t81;
t78 = sin(pkin(12));
t105 = t57 * t78;
t32 = mrSges(6,1) * t105 + mrSges(6,2) * t104;
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t58 = t78 * t87 + t81 * t84;
t23 = t58 * t57;
t56 = -t78 * t84 + t81 * t87;
t24 = t56 * t57;
t9 = t23 * mrSges(7,1) + t24 * mrSges(7,2);
t114 = t32 + t9;
t80 = sin(pkin(6));
t86 = sin(qJ(2));
t103 = t80 * t86;
t83 = cos(pkin(6));
t47 = -t103 * t85 + t83 * t88;
t48 = t103 * t88 + t83 * t85;
t25 = -t82 * t47 + t48 * t79;
t22 = t25 ^ 2;
t101 = -qJ(4) - pkin(8);
t65 = t101 * t88;
t95 = t101 * t85;
t41 = -t65 * t79 - t82 * t95;
t113 = t41 ^ 2;
t77 = t88 ^ 2;
t112 = 0.2e1 * t41;
t110 = t81 / 0.2e1;
t67 = pkin(3) * t79 + qJ(5);
t109 = pkin(9) + t67;
t108 = Ifges(6,4) * t78;
t107 = Ifges(6,4) * t81;
t106 = t25 * t41;
t89 = cos(qJ(2));
t102 = t80 * t89;
t55 = t79 * t85 - t82 * t88;
t71 = -pkin(3) * t88 - pkin(2);
t33 = pkin(4) * t55 - qJ(5) * t57 + t71;
t43 = -t82 * t65 + t79 * t95;
t11 = t78 * t33 + t81 * t43;
t100 = Ifges(7,5) * t58 + Ifges(7,6) * t56;
t61 = -t81 * mrSges(6,1) + t78 * mrSges(6,2);
t99 = t61 - mrSges(5,1);
t98 = t78 ^ 2 + t81 ^ 2;
t97 = t85 ^ 2 + t77;
t96 = Ifges(7,5) * t24 - Ifges(7,6) * t23 + Ifges(7,3) * t55;
t70 = -pkin(3) * t82 - pkin(4);
t36 = t55 * mrSges(5,1) + t57 * mrSges(5,2);
t10 = t81 * t33 - t43 * t78;
t37 = -t56 * mrSges(7,1) + t58 * mrSges(7,2);
t94 = -t10 * t78 + t11 * t81;
t27 = t47 * t79 + t48 * t82;
t14 = -t102 * t81 - t27 * t78;
t15 = -t102 * t78 + t27 * t81;
t93 = -t14 * t78 + t15 * t81;
t92 = -t47 * t85 + t48 * t88;
t74 = t80 ^ 2;
t68 = t74 * t89 ^ 2;
t64 = -mrSges(4,1) * t88 + mrSges(4,2) * t85;
t63 = Ifges(6,1) * t78 + t107;
t62 = Ifges(6,2) * t81 + t108;
t60 = -pkin(5) * t81 + t70;
t50 = t109 * t81;
t49 = t109 * t78;
t39 = Ifges(7,1) * t58 + Ifges(7,4) * t56;
t38 = Ifges(7,4) * t58 + Ifges(7,2) * t56;
t35 = mrSges(6,1) * t55 - mrSges(6,3) * t104;
t34 = -mrSges(6,2) * t55 - mrSges(6,3) * t105;
t31 = -t49 * t84 + t50 * t87;
t30 = -t49 * t87 - t50 * t84;
t18 = pkin(5) * t105 + t41;
t17 = Ifges(6,5) * t55 + (Ifges(6,1) * t81 - t108) * t57;
t16 = Ifges(6,6) * t55 + (-Ifges(6,2) * t78 + t107) * t57;
t13 = mrSges(7,1) * t55 - mrSges(7,3) * t24;
t12 = -mrSges(7,2) * t55 - mrSges(7,3) * t23;
t8 = -pkin(9) * t105 + t11;
t7 = Ifges(7,1) * t24 - Ifges(7,4) * t23 + Ifges(7,5) * t55;
t6 = Ifges(7,4) * t24 - Ifges(7,2) * t23 + Ifges(7,6) * t55;
t5 = pkin(5) * t55 - pkin(9) * t104 + t10;
t4 = t14 * t84 + t15 * t87;
t3 = t14 * t87 - t15 * t84;
t2 = t5 * t84 + t8 * t87;
t1 = t5 * t87 - t8 * t84;
t19 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t22) + m(5) * (t27 ^ 2 + t22 + t68) + m(6) * (t14 ^ 2 + t15 ^ 2 + t22) + m(4) * (t47 ^ 2 + t48 ^ 2 + t68) + m(3) * (t74 * t86 ^ 2 + t83 ^ 2 + t68); -t27 * t55 * mrSges(5,3) + t4 * t12 + t3 * t13 + t14 * t35 + t15 * t34 + t92 * mrSges(4,3) + (t57 * mrSges(5,3) + t114) * t25 + (-t86 * mrSges(3,2) + (mrSges(3,1) - t36 - t64) * t89) * t80 + m(7) * (t1 * t3 + t18 * t25 + t2 * t4) + m(5) * (-t102 * t71 + t27 * t43 + t106) + m(6) * (t10 * t14 + t11 * t15 + t106) + m(4) * (pkin(2) * t102 + pkin(8) * t92); Ifges(4,2) * t77 - 0.2e1 * pkin(2) * t64 + 0.2e1 * t1 * t13 + 0.2e1 * t10 * t35 + 0.2e1 * t11 * t34 + 0.2e1 * t2 * t12 + 0.2e1 * t18 * t9 - t23 * t6 + t24 * t7 + t32 * t112 + 0.2e1 * t71 * t36 + Ifges(3,3) + (Ifges(4,1) * t85 + 0.2e1 * Ifges(4,4) * t88) * t85 + 0.2e1 * t97 * pkin(8) * mrSges(4,3) + (mrSges(5,3) * t112 + Ifges(5,1) * t57 - t78 * t16 + t81 * t17) * t57 + m(4) * (pkin(8) ^ 2 * t97 + pkin(2) ^ 2) + m(5) * (t43 ^ 2 + t71 ^ 2 + t113) + m(6) * (t10 ^ 2 + t11 ^ 2 + t113) + m(7) * (t1 ^ 2 + t18 ^ 2 + t2 ^ 2) + (-0.2e1 * t43 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t55 + (Ifges(6,5) * t81 - Ifges(6,6) * t78 - (2 * Ifges(5,4))) * t57 + t96) * t55; t47 * mrSges(4,1) - t48 * mrSges(4,2) - t27 * mrSges(5,2) + (-t3 * t58 + t4 * t56) * mrSges(7,3) + t93 * mrSges(6,3) + (t37 + t99) * t25 + m(7) * (t25 * t60 + t3 * t30 + t31 * t4) + m(6) * (t25 * t70 + t67 * t93) + m(5) * (-t25 * t82 + t27 * t79) * pkin(3); m(7) * (t1 * t30 + t18 * t60 + t2 * t31) + (-t85 * mrSges(4,1) - t88 * mrSges(4,2)) * pkin(8) + (-t1 * t58 + t2 * t56) * mrSges(7,3) + (m(5) * (-t41 * t82 + t43 * t79) + (-t79 * t55 - t82 * t57) * mrSges(5,3)) * pkin(3) + (t63 * t110 - t78 * t62 / 0.2e1 + Ifges(5,5)) * t57 + t99 * t41 + t94 * mrSges(6,3) + m(6) * (t41 * t70 + t67 * t94) + Ifges(4,5) * t85 + Ifges(4,6) * t88 + t78 * t17 / 0.2e1 + t70 * t32 + t56 * t6 / 0.2e1 + t58 * t7 / 0.2e1 + t60 * t9 - Ifges(5,6) * t55 + t18 * t37 - t23 * t38 / 0.2e1 + t24 * t39 / 0.2e1 - t43 * mrSges(5,2) + t30 * t13 + t31 * t12 + t16 * t110 + (Ifges(6,5) * t78 + Ifges(6,6) * t81 + t100) * t55 / 0.2e1 + (t81 * t34 - t78 * t35) * t67; 0.2e1 * t60 * t37 + t56 * t38 + t58 * t39 + 0.2e1 * t70 * t61 + t81 * t62 + t78 * t63 + Ifges(4,3) + Ifges(5,3) + m(7) * (t30 ^ 2 + t31 ^ 2 + t60 ^ 2) + m(6) * (t67 ^ 2 * t98 + t70 ^ 2) + m(5) * (t79 ^ 2 + t82 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t82 - mrSges(5,2) * t79) * pkin(3) + 0.2e1 * (-t30 * t58 + t31 * t56) * mrSges(7,3) + 0.2e1 * t98 * t67 * mrSges(6,3); -m(5) * t102 + m(7) * (t3 * t56 + t4 * t58) + m(6) * (t14 * t81 + t15 * t78); t58 * t12 + t56 * t13 + t78 * t34 + t81 * t35 + m(7) * (t1 * t56 + t2 * t58) + m(6) * (t10 * t81 + t11 * t78) + m(5) * t71 + t36; m(7) * (t30 * t56 + t31 * t58); m(5) + m(6) * t98 + m(7) * (t56 ^ 2 + t58 ^ 2); 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t25; m(6) * t41 + m(7) * t18 + t114; m(6) * t70 + m(7) * t60 + t37 + t61; 0; m(6) + m(7); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t96; mrSges(7,1) * t30 - t31 * mrSges(7,2) + t100; -t37; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
