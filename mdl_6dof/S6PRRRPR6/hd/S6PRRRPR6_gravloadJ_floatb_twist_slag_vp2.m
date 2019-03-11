% Calculate Gravitation load on the joints for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:13
% EndTime: 2019-03-08 23:31:15
% DurationCPUTime: 1.06s
% Computational Cost: add. (692->126), mult. (1784->183), div. (0->0), fcn. (2182->12), ass. (0->76)
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t123 = pkin(4) * t64 + qJ(5) * t60;
t59 = sin(qJ(6));
t63 = cos(qJ(6));
t109 = -t59 * mrSges(7,1) - t63 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t108 = -m(7) * pkin(5) - mrSges(7,1) * t63 + mrSges(7,2) * t59 - mrSges(5,1) - mrSges(6,1);
t122 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t121 = m(6) + m(7);
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t120 = t65 * mrSges(4,1) - mrSges(4,2) * t61 + mrSges(3,1);
t119 = mrSges(3,2) - mrSges(4,3);
t57 = sin(pkin(11));
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t92 = cos(pkin(11));
t93 = cos(pkin(6));
t74 = t93 * t92;
t43 = t57 * t66 + t62 * t74;
t58 = sin(pkin(6));
t82 = t58 * t92;
t20 = -t43 * t61 - t65 * t82;
t117 = t123 * t20;
t100 = t58 * t65;
t83 = t57 * t93;
t45 = -t62 * t83 + t66 * t92;
t22 = t100 * t57 - t45 * t61;
t116 = t123 * t22;
t101 = t58 * t62;
t46 = -t101 * t61 + t65 * t93;
t115 = t123 * t46;
t114 = pkin(4) * t121 - t108;
t113 = t108 * t64 + t109 * t60 - mrSges(4,1);
t112 = mrSges(4,2) - m(7) * (pkin(9) - pkin(10)) + t122;
t111 = m(7) * pkin(10) + t122;
t105 = pkin(3) * t65;
t42 = t57 * t62 - t66 * t74;
t103 = t42 * t61;
t44 = t62 * t92 + t66 * t83;
t102 = t44 * t61;
t99 = t58 * t66;
t98 = t60 * t65;
t97 = t64 * t65;
t96 = t65 * t66;
t95 = pkin(2) * t99 + pkin(8) * t101;
t91 = t61 * t99;
t90 = t60 * t99;
t88 = -t42 * pkin(2) + pkin(8) * t43;
t87 = -t44 * pkin(2) + pkin(8) * t45;
t17 = t20 * pkin(3);
t21 = t43 * t65 - t61 * t82;
t86 = pkin(9) * t21 + t17;
t18 = t22 * pkin(3);
t23 = t57 * t58 * t61 + t45 * t65;
t85 = pkin(9) * t23 + t18;
t41 = t46 * pkin(3);
t47 = t100 * t62 + t61 * t93;
t84 = pkin(9) * t47 + t41;
t81 = t58 * pkin(3) * t96 + pkin(9) * t91 + t95;
t76 = -pkin(9) * t103 - t42 * t105 + t88;
t75 = -pkin(9) * t102 - t44 * t105 + t87;
t68 = -qJ(5) * t121 + t109;
t28 = (t60 * t62 + t64 * t96) * t58;
t27 = -t101 * t64 + t65 * t90;
t25 = t47 * t64 - t90;
t24 = t47 * t60 + t64 * t99;
t12 = -t44 * t97 + t45 * t60;
t11 = -t44 * t98 - t45 * t64;
t10 = -t42 * t97 + t43 * t60;
t9 = -t42 * t98 - t43 * t64;
t6 = t23 * t64 + t44 * t60;
t5 = t23 * t60 - t44 * t64;
t4 = t21 * t64 + t42 * t60;
t3 = t21 * t60 - t42 * t64;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t121) * g(3) (-m(4) * t95 - m(5) * t81 - t121 * (t28 * pkin(4) + t27 * qJ(5) + t81) + (t119 * t62 - t120 * t66) * t58 + t108 * t28 + t109 * t27 + t111 * t91) * g(3) + (-m(4) * t88 - m(5) * t76 - t121 * (t10 * pkin(4) + qJ(5) * t9 + t76) + t119 * t43 + t120 * t42 + t109 * t9 + t108 * t10 - t111 * t103) * g(2) + (-m(4) * t87 - m(5) * t75 - t121 * (t12 * pkin(4) + qJ(5) * t11 + t75) + t119 * t45 + t120 * t44 + t108 * t12 + t109 * t11 - t111 * t102) * g(1) (-m(5) * t84 - m(6) * (t84 + t115) - m(7) * (t41 + t115) + t112 * t47 + t113 * t46) * g(3) + (-m(5) * t86 - m(6) * (t86 + t117) - m(7) * (t17 + t117) + t112 * t21 + t113 * t20) * g(2) + (-m(5) * t85 - m(6) * (t85 + t116) - m(7) * (t18 + t116) + t112 * t23 + t113 * t22) * g(1) (t114 * t24 + t25 * t68) * g(3) + (t114 * t3 + t4 * t68) * g(2) + (t114 * t5 + t6 * t68) * g(1), t121 * (-g(1) * t5 - g(2) * t3 - g(3) * t24) -g(1) * ((t5 * t63 - t59 * t6) * mrSges(7,1) + (-t5 * t59 - t6 * t63) * mrSges(7,2)) - g(2) * ((t3 * t63 - t4 * t59) * mrSges(7,1) + (-t3 * t59 - t4 * t63) * mrSges(7,2)) - g(3) * ((t24 * t63 - t25 * t59) * mrSges(7,1) + (-t24 * t59 - t25 * t63) * mrSges(7,2))];
taug  = t1(:);
