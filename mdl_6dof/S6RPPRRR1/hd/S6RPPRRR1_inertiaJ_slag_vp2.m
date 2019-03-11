% Calculate joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:58
% EndTime: 2019-03-09 02:18:00
% DurationCPUTime: 0.72s
% Computational Cost: add. (1468->169), mult. (2682->247), div. (0->0), fcn. (2964->10), ass. (0->77)
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t90 = t67 ^ 2 + t70 ^ 2;
t116 = mrSges(7,3) * t90;
t63 = sin(pkin(11));
t65 = cos(pkin(11));
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t42 = -t69 * t63 + t72 * t65;
t43 = t72 * t63 + t69 * t65;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t34 = -t71 * t42 + t68 * t43;
t36 = t68 * t42 + t71 * t43;
t97 = t67 * mrSges(7,3);
t15 = -t34 * mrSges(7,2) - t36 * t97;
t98 = t36 * t70;
t16 = t34 * mrSges(7,1) - mrSges(7,3) * t98;
t115 = t70 * t15 - t67 * t16;
t47 = -t70 * mrSges(7,1) + t67 * mrSges(7,2);
t114 = t36 * t116 + t34 * t47;
t81 = mrSges(7,1) * t67 + mrSges(7,2) * t70;
t14 = t81 * t36;
t64 = sin(pkin(10));
t51 = t64 * pkin(1) + qJ(3);
t103 = pkin(7) + t51;
t87 = t103 * t63;
t88 = t65 * t103;
t27 = -t69 * t87 + t72 * t88;
t18 = t42 * pkin(8) + t27;
t26 = -t69 * t88 - t72 * t87;
t77 = -t43 * pkin(8) + t26;
t8 = t68 * t18 - t71 * t77;
t113 = m(7) * t8 + t14;
t10 = t71 * t18 + t68 * t77;
t106 = pkin(5) * t34;
t66 = cos(pkin(10));
t53 = -t66 * pkin(1) - pkin(2);
t46 = -t65 * pkin(3) + t53;
t37 = -t42 * pkin(4) + t46;
t13 = -t36 * pkin(9) + t106 + t37;
t3 = t70 * t10 + t67 * t13;
t105 = t3 * t70;
t2 = -t67 * t10 + t70 * t13;
t82 = -t2 * t67 + t105;
t112 = m(7) * t82 + t115;
t111 = t8 ^ 2;
t110 = 0.2e1 * t8;
t109 = t34 ^ 2;
t108 = t42 ^ 2;
t60 = t65 ^ 2;
t107 = 0.2e1 * t42;
t104 = t8 * t34;
t101 = Ifges(7,4) * t67;
t100 = Ifges(7,4) * t70;
t99 = t36 * t67;
t94 = Ifges(7,5) * t98 + Ifges(7,3) * t34;
t93 = t34 * mrSges(6,1) + t36 * mrSges(6,2);
t92 = Ifges(7,5) * t67 + Ifges(7,6) * t70;
t91 = t63 ^ 2 + t60;
t48 = Ifges(7,2) * t70 + t101;
t49 = Ifges(7,1) * t67 + t100;
t89 = t70 * t48 + t67 * t49 + Ifges(6,3);
t86 = t90 * pkin(9);
t54 = t68 * pkin(4) + pkin(9);
t85 = t90 * t54;
t84 = -t65 * mrSges(4,1) + t63 * mrSges(4,2);
t83 = -t42 * mrSges(5,1) + t43 * mrSges(5,2);
t80 = 0.2e1 * t116;
t79 = -t83 - t93;
t78 = (t71 * mrSges(6,1) - t68 * mrSges(6,2)) * pkin(4);
t11 = Ifges(7,6) * t34 + (-Ifges(7,2) * t67 + t100) * t36;
t12 = Ifges(7,5) * t34 + (Ifges(7,1) * t70 - t101) * t36;
t76 = -t10 * mrSges(6,2) + mrSges(7,3) * t105 - t2 * t97 - t48 * t99 / 0.2e1 + t49 * t98 / 0.2e1 + Ifges(6,5) * t36 + t67 * t12 / 0.2e1 + t70 * t11 / 0.2e1 + (t47 - mrSges(6,1)) * t8 + (t92 / 0.2e1 - Ifges(6,6)) * t34;
t55 = -t71 * pkin(4) - pkin(5);
t33 = t36 ^ 2;
t1 = [Ifges(4,2) * t60 + 0.2e1 * t53 * t84 + Ifges(5,2) * t108 + 0.2e1 * t46 * t83 + 0.2e1 * t37 * t93 + t14 * t110 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + t27 * mrSges(5,3) * t107 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t63 + 0.2e1 * Ifges(4,4) * t65) * t63 + (-0.2e1 * t26 * mrSges(5,3) + Ifges(5,1) * t43 + Ifges(5,4) * t107) * t43 + (-0.2e1 * t10 * mrSges(6,3) + Ifges(6,2) * t34 + t94) * t34 + (mrSges(6,3) * t110 + Ifges(6,1) * t36 - t67 * t11 + t70 * t12 + (-Ifges(7,6) * t67 - (2 * Ifges(6,4))) * t34) * t36 + m(4) * (t91 * t51 ^ 2 + t53 ^ 2) + m(6) * (t10 ^ 2 + t37 ^ 2 + t111) + m(5) * (t26 ^ 2 + t27 ^ 2 + t46 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t111) + m(3) * (t64 ^ 2 + t66 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t66 * mrSges(3,1) - t64 * mrSges(3,2)) * pkin(1) + 0.2e1 * t91 * t51 * mrSges(4,3); t34 * t14 + t115 * t36 + m(7) * (t82 * t36 + t104) + m(6) * (t10 * t36 + t104) + m(5) * (t26 * t42 + t27 * t43); m(3) + m(7) * (t90 * t33 + t109) + m(6) * (t33 + t109) + m(5) * (t43 ^ 2 + t108) + m(4) * t91; t67 * t15 + t70 * t16 + m(7) * (t70 * t2 + t67 * t3) + m(6) * t37 + m(5) * t46 + m(4) * t53 - t79 + t84; 0; m(7) * t90 + m(4) + m(5) + m(6); (m(6) * (t10 * t68 - t71 * t8) + (-t68 * t34 - t71 * t36) * mrSges(6,3)) * pkin(4) + t76 + Ifges(5,6) * t42 + Ifges(5,5) * t43 + t26 * mrSges(5,1) - t27 * mrSges(5,2) + t113 * t55 + t112 * t54; m(7) * (t55 * t34 + t36 * t85) + m(6) * (-t34 * t71 + t36 * t68) * pkin(4) + t79 + t114; 0; 0.2e1 * t55 * t47 + Ifges(5,3) + 0.2e1 * t78 + t54 * t80 + m(7) * (t90 * t54 ^ 2 + t55 ^ 2) + m(6) * (t68 ^ 2 + t71 ^ 2) * pkin(4) ^ 2 + t89; -t113 * pkin(5) + t112 * pkin(9) + t76; m(7) * (t36 * t86 - t106) - t93 + t114; 0; m(7) * (-pkin(5) * t55 + pkin(9) * t85) + (-pkin(5) + t55) * t47 + t78 + (t85 + t86) * mrSges(7,3) + t89; -0.2e1 * pkin(5) * t47 + m(7) * (t90 * pkin(9) ^ 2 + pkin(5) ^ 2) + pkin(9) * t80 + t89; t2 * mrSges(7,1) - t3 * mrSges(7,2) - Ifges(7,6) * t99 + t94; -t14; -t47; -t81 * t54 + t92; -t81 * pkin(9) + t92; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
