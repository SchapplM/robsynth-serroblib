% Calculate Gravitation load on the joints for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:26:38
% EndTime: 2018-11-23 17:26:39
% DurationCPUTime: 1.12s
% Computational Cost: add. (1492->153), mult. (1402->194), div. (0->0), fcn. (1354->18), ass. (0->81)
t153 = mrSges(6,2) - mrSges(7,3);
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t152 = t81 * mrSges(7,1) - t78 * mrSges(7,2) + mrSges(6,1);
t154 = m(7) * pkin(5) + t152;
t106 = -m(7) * pkin(11) + t153;
t74 = sin(pkin(6));
t80 = sin(qJ(1));
t132 = t74 * t80;
t136 = cos(qJ(1));
t82 = cos(qJ(2));
t117 = t136 * t82;
t123 = pkin(6) - qJ(2);
t108 = sin(t123);
t122 = pkin(6) + qJ(2);
t107 = sin(t122);
t99 = t107 / 0.2e1;
t92 = t99 - t108 / 0.2e1;
t41 = -t80 * t92 + t117;
t72 = pkin(12) + qJ(4);
t66 = sin(t72);
t67 = cos(t72);
t19 = t67 * t132 - t41 * t66;
t101 = cos(t122) / 0.2e1;
t109 = cos(t123);
t54 = t101 - t109 / 0.2e1;
t76 = cos(pkin(6));
t150 = t54 * t66 + t67 * t76;
t75 = cos(pkin(12));
t65 = t75 * pkin(3) + pkin(2);
t73 = sin(pkin(12));
t149 = -m(4) * pkin(2) - m(5) * t65 - mrSges(4,1) * t75 + mrSges(4,2) * t73 - mrSges(3,1);
t148 = -m(4) - m(5);
t147 = m(6) + m(7);
t77 = -pkin(9) - qJ(3);
t86 = -m(4) * qJ(3) + m(5) * t77 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t145 = t78 * mrSges(7,1) + t81 * mrSges(7,2) - t86;
t68 = qJ(5) + t72;
t63 = sin(t68);
t64 = cos(t68);
t30 = t54 * t63 + t64 * t76;
t31 = -t54 * t64 + t63 * t76;
t144 = -t152 * t30 + t153 * t31;
t17 = -t64 * t132 + t41 * t63;
t18 = t63 * t132 + t41 * t64;
t143 = t152 * t17 + t153 * t18;
t118 = t74 * t136;
t130 = t80 * t82;
t38 = t136 * t92 + t130;
t13 = -t64 * t118 - t38 * t63;
t14 = -t63 * t118 + t38 * t64;
t142 = -t152 * t13 + t153 * t14;
t141 = t67 * mrSges(5,1) - t66 * mrSges(5,2) - t106 * t63 + t154 * t64 - t149;
t124 = t136 * pkin(1) + pkin(8) * t132;
t115 = -t80 * pkin(1) + pkin(8) * t118;
t114 = t13 * pkin(5) + t14 * pkin(11);
t113 = -t17 * pkin(5) + pkin(11) * t18;
t112 = t30 * pkin(5) + pkin(11) * t31;
t110 = t66 * t118 - t38 * t67;
t105 = t73 * t118;
t104 = t19 * pkin(4);
t103 = t150 * pkin(4);
t79 = sin(qJ(2));
t84 = t109 / 0.2e1 + t101;
t40 = t136 * t79 + t80 * t84;
t48 = pkin(4) * t67 + t65;
t55 = pkin(3) * t73 + pkin(4) * t66;
t71 = -pkin(10) + t77;
t102 = t55 * t132 - t40 * t71 + t41 * t48 + t124;
t100 = t108 / 0.2e1;
t93 = t100 - t107 / 0.2e1;
t88 = t67 * t118 + t38 * t66;
t85 = t88 * pkin(4);
t53 = t99 + t100;
t42 = t80 * t93 + t117;
t39 = -t136 * t93 + t130;
t37 = -t136 * t84 + t79 * t80;
t20 = t66 * t132 + t41 * t67;
t2 = t18 * t81 + t40 * t78;
t1 = -t18 * t78 + t40 * t81;
t3 = [(t80 * mrSges(2,1) + t136 * mrSges(2,2) - m(3) * t115 + t38 * mrSges(3,1) - mrSges(3,3) * t118 - m(4) * (-pkin(2) * t38 + t115) - (-t38 * t75 + t105) * mrSges(4,1) - (t75 * t118 + t38 * t73) * mrSges(4,2) - m(5) * (pkin(3) * t105 - t38 * t65 + t115) - t110 * mrSges(5,1) - t88 * mrSges(5,2) + t106 * t13 + t154 * t14 + t145 * t37 + t147 * (-t55 * t118 - t37 * t71 + t38 * t48 - t115)) * g(1) + (-t136 * mrSges(2,1) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t102 - t18 * mrSges(6,1) - m(7) * (t18 * pkin(5) + t102) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t149 * t41 + (mrSges(2,2) + (-mrSges(4,2) * t75 - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t73) * t74) * t80 + t106 * t17 + t86 * t40 + (-m(3) + t148) * t124) * g(2) (-t147 * (t53 * t48 + t54 * t71) + t145 * t54 - t141 * t53) * g(3) + (-t147 * (-t37 * t48 - t39 * t71) - t145 * t39 + t141 * t37) * g(2) + (-t147 * (-t40 * t48 - t42 * t71) - t145 * t42 + t141 * t40) * g(1) (-g(1) * t40 - g(2) * t37 + g(3) * t53) * (t147 - t148) (-t150 * mrSges(5,1) - (t54 * t67 - t66 * t76) * mrSges(5,2) - m(6) * t103 - m(7) * (t103 + t112) + t144) * g(3) + (t88 * mrSges(5,1) - t110 * mrSges(5,2) + m(6) * t85 - m(7) * (t114 - t85) + t142) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t104 - m(7) * (t104 + t113) + t143) * g(1) (-m(7) * t112 + t144) * g(3) + (-m(7) * t114 + t142) * g(2) + (-m(7) * t113 + t143) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t78 + t37 * t81) * mrSges(7,1) + (-t14 * t81 - t37 * t78) * mrSges(7,2)) - g(3) * ((-t31 * t78 - t53 * t81) * mrSges(7,1) + (-t31 * t81 + t53 * t78) * mrSges(7,2))];
taug  = t3(:);
