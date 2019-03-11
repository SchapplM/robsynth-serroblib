% Calculate Gravitation load on the joints for
% S6RRPRRR4
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:50
% EndTime: 2019-03-09 13:27:52
% DurationCPUTime: 1.24s
% Computational Cost: add. (941->149), mult. (1930->203), div. (0->0), fcn. (2365->14), ass. (0->81)
t158 = mrSges(6,2) - mrSges(7,3);
t78 = sin(qJ(6));
t82 = cos(qJ(6));
t157 = mrSges(7,1) * t82 - mrSges(7,2) * t78 + mrSges(6,1);
t159 = -m(7) * pkin(5) - t157;
t107 = -m(7) * pkin(11) + t158;
t152 = m(6) + m(7);
t153 = -m(4) - m(5);
t117 = t152 - t153;
t124 = cos(pkin(12));
t75 = sin(pkin(12));
t80 = sin(qJ(2));
t84 = cos(qJ(2));
t54 = -t80 * t124 - t84 * t75;
t85 = cos(qJ(1));
t128 = t84 * t85;
t81 = sin(qJ(1));
t133 = t80 * t81;
t77 = cos(pkin(6));
t156 = t77 * t128 - t133;
t76 = sin(pkin(6));
t137 = t76 * t81;
t125 = t54 * t77;
t94 = t84 * t124 - t80 * t75;
t34 = t125 * t81 + t85 * t94;
t79 = sin(qJ(4));
t83 = cos(qJ(4));
t19 = t83 * t137 - t34 * t79;
t47 = t54 * t76;
t155 = t47 * t79 + t77 * t83;
t151 = m(5) * pkin(3) + mrSges(5,1) * t83 - mrSges(5,2) * t79 + mrSges(4,1);
t99 = m(5) * pkin(9) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t150 = -mrSges(7,1) * t78 - mrSges(7,2) * t82 - t99;
t74 = qJ(4) + qJ(5);
t72 = sin(t74);
t73 = cos(t74);
t40 = t47 * t72 + t73 * t77;
t41 = -t47 * t73 + t72 * t77;
t149 = -t157 * t40 + t158 * t41;
t17 = -t137 * t73 + t34 * t72;
t18 = t137 * t72 + t34 * t73;
t148 = t157 * t17 + t158 * t18;
t136 = t76 * t85;
t29 = t125 * t85 - t81 * t94;
t13 = -t73 * t136 + t29 * t72;
t14 = -t72 * t136 - t29 * t73;
t147 = -t157 * t13 + t158 * t14;
t130 = t81 * t84;
t132 = t80 * t85;
t51 = -t77 * t130 - t132;
t146 = t107 * t72 + t159 * t73 - t151;
t143 = pkin(2) * t84;
t67 = t76 * t143;
t123 = t79 * t137;
t66 = t79 * t136;
t116 = -m(3) * pkin(1) - mrSges(2,1);
t114 = t13 * pkin(5) + pkin(11) * t14;
t113 = -t17 * pkin(5) + pkin(11) * t18;
t112 = t40 * pkin(5) + pkin(11) * t41;
t89 = t77 * t94;
t33 = t54 * t85 - t81 * t89;
t71 = pkin(1) + t143;
t65 = t85 * t71;
t70 = pkin(4) * t83 + pkin(3);
t86 = -pkin(10) - pkin(9);
t109 = pkin(4) * t123 + t33 * t86 + t34 * t70 + t65;
t108 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t106 = t156 * pkin(2);
t105 = t19 * pkin(4);
t104 = t155 * pkin(4);
t97 = -t136 * t83 + t29 * t79;
t95 = t117 * (pkin(2) * t77 * t80 + (-pkin(8) - qJ(3)) * t76) + mrSges(2,2);
t91 = t97 * pkin(4);
t52 = -t133 * t77 + t128;
t50 = -t132 * t77 - t130;
t46 = t94 * t76;
t30 = t81 * t54 + t85 * t89;
t20 = t34 * t83 + t123;
t2 = t18 * t82 - t33 * t78;
t1 = -t18 * t78 - t33 * t82;
t3 = [(-t52 * mrSges(3,1) - t51 * mrSges(3,2) - m(4) * t65 - t34 * mrSges(4,1) - m(5) * (pkin(3) * t34 + t65) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t109 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t109) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t116 * t85 + t107 * t17 + (t108 * t76 + t95) * t81 + t99 * t33) * g(2) + (-t50 * mrSges(3,1) + t156 * mrSges(3,2) - t66 * mrSges(5,1) + t107 * t13 - t151 * t29 + (t117 * t71 - t116) * t81 + ((-mrSges(5,2) * t83 + t108) * t76 + t95) * t85 - t159 * t14 + t150 * t30 + t152 * (-pkin(4) * t66 - t29 * t70 + t30 * t86)) * g(1) (-(mrSges(3,1) * t84 - mrSges(3,2) * t80) * t76 + t153 * t67 - t152 * (t46 * t70 + t47 * t86 + t67) - t150 * t47 + t146 * t46) * g(3) + (-mrSges(3,1) * t156 - mrSges(3,2) * t50 - t152 * (t29 * t86 + t30 * t70 + t106) + t153 * t106 + t146 * t30 - t150 * t29) * g(2) + (mrSges(3,2) * t52 - t152 * (t33 * t70 - t34 * t86) + t146 * t33 + t150 * t34 + (-t117 * pkin(2) - mrSges(3,1)) * t51) * g(1) (-t77 * g(3) + (-g(1) * t81 + g(2) * t85) * t76) * t117 (-t155 * mrSges(5,1) - (t47 * t83 - t77 * t79) * mrSges(5,2) - m(6) * t104 - m(7) * (t104 + t112) + t149) * g(3) + (-t97 * mrSges(5,1) - (t29 * t83 + t66) * mrSges(5,2) - m(6) * t91 - m(7) * (t114 + t91) + t147) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t105 - m(7) * (t105 + t113) + t148) * g(1) (-m(7) * t112 + t149) * g(3) + (-m(7) * t114 + t147) * g(2) + (-m(7) * t113 + t148) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t78 - t30 * t82) * mrSges(7,1) + (-t14 * t82 + t30 * t78) * mrSges(7,2)) - g(3) * ((-t41 * t78 - t46 * t82) * mrSges(7,1) + (-t41 * t82 + t46 * t78) * mrSges(7,2))];
taug  = t3(:);
