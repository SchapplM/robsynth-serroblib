% Calculate Gravitation load on the joints for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:20:22
% EndTime: 2018-11-23 17:20:23
% DurationCPUTime: 0.97s
% Computational Cost: add. (1431->131), mult. (1749->175), div. (0->0), fcn. (1735->14), ass. (0->74)
t139 = m(6) + m(7);
t136 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t100 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t95 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t120 = -mrSges(3,2) + mrSges(4,3);
t138 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t137 = -t100 * t72 + t95 * t68 - mrSges(5,1);
t118 = m(5) + t139;
t135 = pkin(9) * t118 + t138;
t71 = sin(qJ(1));
t74 = cos(qJ(2));
t122 = t71 * t74;
t75 = cos(qJ(1));
t117 = pkin(6) - qJ(2);
t102 = sin(t117);
t116 = pkin(6) + qJ(2);
t101 = sin(t116);
t91 = t101 / 0.2e1;
t83 = t91 - t102 / 0.2e1;
t43 = t75 * t83 + t122;
t66 = sin(pkin(6));
t125 = t66 * t75;
t70 = sin(qJ(2));
t103 = cos(t116);
t93 = cos(t117) / 0.2e1;
t85 = t93 + t103 / 0.2e1;
t42 = t70 * t71 - t75 * t85;
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t82 = t125 * t73 - t42 * t69;
t134 = t43 * t72 + t68 * t82;
t133 = -t43 * t68 + t72 * t82;
t113 = m(4) + t118;
t132 = qJ(3) * t113 + t120;
t131 = -t69 * mrSges(5,1) - t120 - t139 * (-pkin(10) * t73 + qJ(3)) + (-m(4) - m(5)) * qJ(3) - t136 * t73;
t129 = pkin(4) * t69;
t126 = t66 * t71;
t124 = t68 * t69;
t123 = t69 * t72;
t121 = t75 * t74;
t119 = t75 * pkin(1) + pkin(8) * t126;
t46 = -t71 * t83 + t121;
t114 = t46 * pkin(2) + t119;
t112 = -pkin(1) * t71 + pkin(8) * t125;
t36 = t42 * pkin(2);
t111 = -pkin(9) * t42 - t36;
t45 = t75 * t70 + t71 * t85;
t38 = t45 * pkin(2);
t110 = -pkin(9) * t45 - t38;
t92 = t102 / 0.2e1;
t55 = t91 + t92;
t54 = t55 * pkin(2);
t109 = pkin(9) * t55 + t54;
t104 = pkin(3) * t126 + t114;
t99 = -t43 * pkin(2) + t112;
t90 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t66;
t86 = pkin(3) * t125 + t99;
t84 = t92 - t101 / 0.2e1;
t19 = t125 * t69 + t42 * t73;
t80 = -t139 * pkin(10) + t136;
t67 = cos(pkin(6));
t56 = t93 - t103 / 0.2e1;
t47 = t71 * t84 + t121;
t44 = -t75 * t84 + t122;
t41 = -t55 * t69 + t67 * t73;
t40 = -t55 * t73 - t67 * t69;
t18 = t126 * t73 + t45 * t69;
t17 = t126 * t69 - t45 * t73;
t11 = t41 * t68 - t56 * t72;
t2 = t18 * t72 + t46 * t68;
t1 = t18 * t68 - t46 * t72;
t3 = [(-t75 * mrSges(2,1) - m(3) * t119 - m(4) * t114 - m(5) * t104 - t18 * mrSges(5,1) - t132 * t45 - t100 * t2 + t95 * t1 + t90 * t71 - t135 * t46 + t80 * t17 - t139 * (t18 * pkin(4) + t104)) * g(2) + (t71 * mrSges(2,1) - m(3) * t112 - m(4) * t99 - m(5) * t86 - t82 * mrSges(5,1) + t132 * t42 - t100 * t133 + t95 * t134 + t90 * t75 + t135 * t43 + t80 * t19 + t139 * (-pkin(4) * t82 - t86)) * g(1) (-m(4) * t54 - m(5) * t109 - t139 * (t56 * t129 + t109) - t100 * (t123 * t56 + t55 * t68) + t95 * (t124 * t56 - t55 * t72) - t138 * t55 + t131 * t56) * g(3) + (m(4) * t36 - m(5) * t111 - t139 * (t44 * t129 + t111) - t100 * (t123 * t44 - t42 * t68) + t95 * (t124 * t44 + t42 * t72) + t138 * t42 + t131 * t44) * g(2) + (m(4) * t38 - m(5) * t110 - t139 * (t47 * t129 + t110) + t95 * (t124 * t47 + t45 * t72) - t100 * (t123 * t47 - t45 * t68) + t138 * t45 + t131 * t47) * g(1) (-g(1) * t45 - g(2) * t42 + g(3) * t55) * t113 (-t139 * (t40 * pkin(4) + pkin(10) * t41) + t136 * t41 + t137 * t40) * g(3) + (-t139 * (t19 * pkin(4) - pkin(10) * t82) - t136 * t82 + t137 * t19) * g(2) + (-t139 * (-t17 * pkin(4) + pkin(10) * t18) + t136 * t18 - t137 * t17) * g(1) (t95 * (t41 * t72 + t56 * t68) + t100 * t11) * g(3) + (-t100 * t134 - t95 * t133) * g(2) + (t1 * t100 + t2 * t95) * g(1) (-g(1) * t1 + g(2) * t134 - g(3) * t11) * m(7)];
taug  = t3(:);
