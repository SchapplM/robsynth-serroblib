% Calculate Gravitation load on the joints for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:37
% EndTime: 2019-03-09 16:03:41
% DurationCPUTime: 1.55s
% Computational Cost: add. (594->146), mult. (1439->184), div. (0->0), fcn. (1677->10), ass. (0->81)
t115 = m(6) + m(7);
t124 = mrSges(5,2) + mrSges(4,3) - mrSges(3,2) - mrSges(6,3);
t127 = t115 * (pkin(9) - qJ(5)) + t124;
t110 = mrSges(4,2) - mrSges(6,1) - mrSges(5,3);
t126 = -m(7) * pkin(5) + t110;
t125 = -m(7) * pkin(10) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2);
t105 = mrSges(7,3) - t125;
t49 = sin(qJ(6));
t120 = t49 * mrSges(7,2);
t106 = t110 - m(7) * (qJ(4) + pkin(5));
t52 = cos(qJ(6));
t63 = t52 * mrSges(7,1) - t120;
t107 = t106 - t63;
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t118 = mrSges(3,1) + t105 * t53 + (t63 - t126) * t50;
t116 = t49 * mrSges(7,1) + t52 * mrSges(7,2);
t117 = t116 - t127;
t100 = sin(qJ(1));
t101 = cos(qJ(2));
t51 = sin(qJ(2));
t102 = cos(qJ(1));
t89 = cos(pkin(6));
t68 = t89 * t102;
t31 = t100 * t51 - t101 * t68;
t90 = qJ(4) * t50;
t97 = t31 * t53;
t112 = -pkin(3) * t97 - t31 * t90;
t67 = t89 * t100;
t33 = t101 * t67 + t102 * t51;
t96 = t33 * t53;
t111 = -pkin(3) * t96 - t33 * t90;
t104 = pkin(9) * t33;
t103 = t31 * pkin(9);
t99 = t31 * t49;
t98 = t31 * t52;
t48 = sin(pkin(6));
t95 = t48 * t51;
t86 = t48 * t101;
t92 = pkin(2) * t86 + pkin(9) * t95;
t85 = t48 * t100;
t91 = t102 * pkin(1) + pkin(8) * t85;
t34 = t101 * t102 - t51 * t67;
t88 = t34 * pkin(2) + t91;
t87 = t48 * t102;
t84 = t50 * t101;
t83 = t52 * t101;
t82 = t53 * t101;
t25 = t31 * pkin(2);
t32 = t100 * t101 + t51 * t68;
t81 = t32 * pkin(9) - t25;
t27 = t33 * pkin(2);
t80 = t34 * pkin(9) - t27;
t12 = t32 * t53 - t50 * t87;
t11 = t32 * t50 + t53 * t87;
t4 = t11 * pkin(3);
t79 = qJ(4) * t12 - t4;
t16 = t34 * t53 + t50 * t85;
t15 = t34 * t50 - t53 * t85;
t8 = t15 * pkin(3);
t78 = qJ(4) * t16 - t8;
t29 = t50 * t95 - t53 * t89;
t24 = t29 * pkin(3);
t30 = t50 * t89 + t53 * t95;
t77 = qJ(4) * t30 - t24;
t76 = t16 * pkin(3) + t88;
t72 = t48 * t82;
t73 = t48 * qJ(4) * t84 + pkin(3) * t72 + t92;
t71 = -pkin(1) * t100 + pkin(8) * t87;
t62 = -t32 * pkin(2) + t71;
t61 = -pkin(3) * t12 + t62;
t60 = qJ(4) * t15 + t76;
t56 = -qJ(4) * t11 + t61;
t23 = t29 * pkin(4);
t9 = t16 * pkin(4);
t7 = t15 * pkin(4);
t5 = t12 * pkin(4);
t3 = t11 * pkin(4);
t2 = t15 * t52 - t33 * t49;
t1 = -t15 * t49 - t33 * t52;
t6 = [(-t102 * mrSges(2,1) + t100 * mrSges(2,2) - m(3) * t91 - t34 * mrSges(3,1) - mrSges(3,3) * t85 - m(4) * (t88 + t104) - m(5) * (t60 + t104) - m(6) * (t60 + t9) - m(7) * (t76 + t9) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t106 * t15 - t127 * t33 - t105 * t16) * g(2) + (t100 * mrSges(2,1) + t102 * mrSges(2,2) - m(3) * t71 + t32 * mrSges(3,1) - mrSges(3,3) * t87 - m(4) * (t62 - t103) - m(5) * (t56 - t103) - m(6) * (-t5 + t56) - m(7) * (-t5 + t61) - t99 * mrSges(7,1) - t98 * mrSges(7,2) - t107 * t11 + t127 * t31 + t105 * t12) * g(1) (-m(4) * t92 - m(5) * t73 - mrSges(7,3) * t72 - t115 * (pkin(4) * t72 + t73) + (-t50 * t83 * mrSges(7,1) - mrSges(3,1) * t101 + t125 * t82 + (t120 + t126) * t84 + (qJ(5) * t115 + t116 - t124) * t51) * t48) * g(3) + (-m(4) * t81 - m(5) * (t81 + t112) - t115 * (-pkin(4) * t97 + t112 - t25) + t117 * t32 + t118 * t31) * g(2) + (-m(4) * t80 - m(5) * (t80 + t111) - t115 * (-pkin(4) * t96 + t111 - t27) + t117 * t34 + t118 * t33) * g(1) (-m(5) * t77 - m(6) * (-t23 + t77) - m(7) * (-t23 - t24) + t107 * t30 + t105 * t29) * g(3) + (-m(5) * t79 - m(6) * (-t3 + t79) - m(7) * (-t3 - t4) + t107 * t12 + t105 * t11) * g(2) + (-m(5) * t78 - m(6) * (-t7 + t78) - m(7) * (-t7 - t8) + t107 * t16 + t105 * t15) * g(1) (m(5) + t115) * (-g(1) * t15 - g(2) * t11 - g(3) * t29) t115 * (g(1) * t33 + g(2) * t31 - g(3) * t86) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t11 * t49 - t98) * mrSges(7,1) + (-t11 * t52 + t99) * mrSges(7,2)) - g(3) * ((-t29 * t49 + t48 * t83) * mrSges(7,1) + (-t29 * t52 - t49 * t86) * mrSges(7,2))];
taug  = t6(:);
