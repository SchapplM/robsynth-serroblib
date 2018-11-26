% Calculate Gravitation load on the joints for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:17:33
% EndTime: 2018-11-23 17:17:34
% DurationCPUTime: 1.34s
% Computational Cost: add. (1639->142), mult. (1738->190), div. (0->0), fcn. (1730->16), ass. (0->67)
t106 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t105 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t121 = sin(pkin(6));
t133 = cos(qJ(1));
t103 = t133 * t121;
t76 = pkin(11) + qJ(4);
t73 = sin(t76);
t74 = cos(t76);
t119 = pkin(6) + qJ(2);
t101 = sin(t119) / 0.2e1;
t120 = pkin(6) - qJ(2);
t107 = sin(t120);
t60 = t101 - t107 / 0.2e1;
t83 = sin(qJ(1));
t85 = cos(qJ(2));
t92 = t133 * t60 + t83 * t85;
t18 = -t73 * t103 + t74 * t92;
t82 = sin(qJ(2));
t102 = cos(t119) / 0.2e1;
t108 = cos(t120);
t87 = t108 / 0.2e1 + t102;
t45 = -t133 * t87 + t82 * t83;
t81 = sin(qJ(5));
t84 = cos(qJ(5));
t1 = t18 * t81 - t45 * t84;
t143 = t18 * t84 + t45 * t81;
t142 = m(6) + m(7);
t77 = sin(pkin(11));
t78 = cos(pkin(11));
t141 = -m(4) * pkin(2) - t78 * mrSges(4,1) + t77 * mrSges(4,2) - mrSges(3,1);
t140 = mrSges(6,3) + mrSges(7,2);
t112 = t83 * t121;
t139 = t105 * t81 - t106 * t84 - mrSges(5,1);
t138 = mrSges(5,2) - t140;
t137 = t74 * mrSges(5,1) - mrSges(5,2) * t73 - t141;
t136 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t134 = pkin(4) * t74;
t132 = t45 * t73;
t48 = t133 * t82 + t83 * t87;
t129 = t48 * t73;
t59 = t101 + t107 / 0.2e1;
t128 = t59 * t73;
t127 = t74 * t81;
t126 = t74 * t84;
t72 = pkin(3) * t78 + pkin(2);
t80 = -pkin(9) - qJ(3);
t125 = -t45 * t72 - t80 * t92;
t91 = t133 * t85 - t83 * t60;
t124 = -t48 * t72 - t80 * t91;
t61 = t102 - t108 / 0.2e1;
t123 = t59 * t72 + t61 * t80;
t122 = t133 * pkin(1) + pkin(8) * t112;
t117 = -pkin(1) * t83 + pkin(8) * t103;
t17 = -t74 * t103 - t73 * t92;
t104 = t77 * pkin(3) * t112 - t48 * t80 + t72 * t91 + t122;
t94 = t77 * t103;
t93 = pkin(3) * t94 + t45 * t80 - t72 * t92 + t117;
t90 = -pkin(10) * t142 + t138;
t79 = cos(pkin(6));
t35 = -t61 * t74 + t73 * t79;
t34 = t61 * t73 + t74 * t79;
t22 = t112 * t73 + t74 * t91;
t21 = -t112 * t74 + t73 * t91;
t11 = t35 * t81 + t59 * t84;
t6 = t22 * t84 + t48 * t81;
t5 = t22 * t81 - t48 * t84;
t2 = [(-t133 * mrSges(2,1) - m(5) * t104 - t22 * mrSges(5,1) + t141 * t91 + t83 * mrSges(2,2) - t106 * t6 + t105 * t5 + t136 * t48 + t90 * t21 - t142 * (t22 * pkin(4) + t104) + (-m(3) - m(4)) * t122 + (-mrSges(4,1) * t77 - t78 * mrSges(4,2) - mrSges(3,3)) * t112) * g(2) + (t83 * mrSges(2,1) + t133 * mrSges(2,2) - m(3) * t117 + t92 * mrSges(3,1) - mrSges(3,3) * t103 - m(4) * (-pkin(2) * t92 + t117) - (-t78 * t92 + t94) * mrSges(4,1) - (t103 * t78 + t77 * t92) * mrSges(4,2) - m(5) * t93 + t18 * mrSges(5,1) + t106 * t143 - t105 * t1 - t136 * t45 + t90 * t17 + t142 * (pkin(4) * t18 - t93)) * g(1) (-m(5) * t123 - t142 * (pkin(10) * t128 + t59 * t134 + t123) - t106 * (t126 * t59 - t61 * t81) + t105 * (t127 * t59 + t61 * t84) - t140 * t128 - t136 * t61 - t137 * t59) * g(3) + (-m(5) * t125 - t142 * (-pkin(10) * t132 - t45 * t134 + t125) - t106 * (-t126 * t45 + t81 * t92) + t105 * (-t127 * t45 - t84 * t92) + t140 * t132 + t136 * t92 + t137 * t45) * g(2) + (-m(5) * t124 - t142 * (-pkin(10) * t129 - t48 * t134 + t124) + t105 * (-t127 * t48 - t84 * t91) + t140 * t129 - t106 * (-t126 * t48 + t81 * t91) + t136 * t91 + t137 * t48) * g(1) (-g(1) * t48 - g(2) * t45 + g(3) * t59) * (m(4) + m(5) + t142) (-t142 * (t34 * pkin(4) + pkin(10) * t35) + t138 * t35 + t139 * t34) * g(3) + (-t142 * (t17 * pkin(4) + pkin(10) * t18) + t138 * t18 + t139 * t17) * g(2) + (-t142 * (-t21 * pkin(4) + pkin(10) * t22) + t138 * t22 - t139 * t21) * g(1) (t105 * (t35 * t84 - t59 * t81) + t106 * t11) * g(3) + (t106 * t1 + t105 * t143) * g(2) + (t105 * t6 + t106 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
