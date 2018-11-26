% Calculate Gravitation load on the joints for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:49:49
% EndTime: 2018-11-23 17:49:50
% DurationCPUTime: 1.26s
% Computational Cost: add. (1628->139), mult. (2019->188), div. (0->0), fcn. (2032->14), ass. (0->77)
t144 = mrSges(4,2) - mrSges(5,3);
t145 = m(6) + m(7);
t149 = m(5) + t145;
t89 = -qJ(4) * t149 + t144;
t142 = mrSges(4,1) - mrSges(5,2) + mrSges(7,2) + mrSges(6,3);
t104 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t148 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t138 = cos(qJ(1));
t72 = sin(pkin(6));
t117 = t72 * t138;
t76 = sin(qJ(1));
t79 = cos(qJ(2));
t129 = t76 * t79;
t122 = pkin(6) - qJ(2);
t106 = sin(t122);
t121 = pkin(6) + qJ(2);
t105 = sin(t121);
t98 = t105 / 0.2e1;
t91 = t98 - t106 / 0.2e1;
t49 = t138 * t91 + t129;
t74 = sin(qJ(3));
t78 = cos(qJ(3));
t20 = t117 * t78 + t49 * t74;
t75 = sin(qJ(2));
t100 = cos(t121) / 0.2e1;
t107 = cos(t122);
t82 = t107 / 0.2e1 + t100;
t48 = -t138 * t82 + t75 * t76;
t73 = sin(qJ(5));
t77 = cos(qJ(5));
t147 = t20 * t73 + t48 * t77;
t146 = t20 * t77 - t48 * t73;
t120 = mrSges(3,2) - mrSges(5,1) - mrSges(4,3);
t143 = pkin(9) * (m(4) + t149) - t120;
t141 = t142 * t78 - t144 * t74 + mrSges(3,1);
t140 = -t104 * t73 + t148 * t77 + t89;
t135 = t48 * t78;
t51 = t138 * t75 + t76 * t82;
t134 = t51 * t78;
t99 = t106 / 0.2e1;
t62 = t98 + t99;
t133 = t62 * t78;
t132 = t72 * t76;
t131 = t73 * t74;
t130 = t74 * t77;
t126 = t138 * pkin(1) + pkin(8) * t132;
t125 = qJ(4) * t74;
t124 = cos(pkin(6));
t116 = t138 * t79;
t52 = -t76 * t91 + t116;
t119 = t52 * pkin(2) + t126;
t115 = -pkin(1) * t76 + pkin(8) * t117;
t81 = t99 - t105 / 0.2e1;
t50 = -t138 * t81 + t129;
t114 = -t48 * pkin(2) + pkin(9) * t50;
t53 = t76 * t81 + t116;
t113 = -t51 * pkin(2) + pkin(9) * t53;
t63 = t100 - t107 / 0.2e1;
t112 = t62 * pkin(2) - pkin(9) * t63;
t21 = -t74 * t117 + t49 * t78;
t25 = t132 * t74 + t52 * t78;
t108 = t25 * pkin(3) + t119;
t103 = -t49 * pkin(2) + t115;
t95 = -pkin(3) * t21 + t103;
t94 = -pkin(3) * t135 - t48 * t125 + t114;
t93 = -pkin(3) * t134 - t51 * t125 + t113;
t92 = pkin(3) * t133 + t62 * t125 + t112;
t84 = -pkin(10) * t145 - t142;
t46 = -t124 * t78 - t63 * t74;
t37 = t46 * pkin(3);
t24 = -t132 * t78 + t52 * t74;
t18 = t24 * pkin(3);
t16 = t20 * pkin(3);
t11 = -t46 * t77 - t62 * t73;
t6 = t24 * t73 + t51 * t77;
t5 = -t24 * t77 + t51 * t73;
t1 = [(-t138 * mrSges(2,1) - m(3) * t126 - t52 * mrSges(3,1) - m(4) * t119 - m(5) * t108 + (-mrSges(3,3) * t72 + mrSges(2,2)) * t76 - t104 * t6 - t148 * t5 + t89 * t24 - t143 * t51 + t84 * t25 - t145 * (t51 * pkin(4) + t108)) * g(2) + (t76 * mrSges(2,1) + t138 * mrSges(2,2) - m(3) * t115 + t49 * mrSges(3,1) - mrSges(3,3) * t117 - m(4) * t103 - m(5) * t95 + t104 * t147 - t148 * t146 - t89 * t20 + t143 * t48 - t84 * t21 + t145 * (t48 * pkin(4) - t95)) * g(1) (-m(4) * t112 - m(5) * t92 - t145 * (-t63 * pkin(4) + pkin(10) * t133 + t92) - t104 * (t131 * t62 - t63 * t77) - t148 * (-t130 * t62 - t63 * t73) - t120 * t63 - t141 * t62) * g(3) + (-m(4) * t114 - m(5) * t94 - t145 * (t50 * pkin(4) - pkin(10) * t135 + t94) - t104 * (-t131 * t48 + t50 * t77) - t148 * (t130 * t48 + t50 * t73) + t120 * t50 + t141 * t48) * g(2) + (-m(4) * t113 - m(5) * t93 - t145 * (t53 * pkin(4) - pkin(10) * t134 + t93) - t148 * (t130 * t51 + t53 * t73) - t104 * (-t131 * t51 + t53 * t77) + t120 * t53 + t141 * t51) * g(1) (m(5) * t37 - t145 * (-pkin(10) * t46 - t37) + t140 * (t124 * t74 - t63 * t78) + t142 * t46) * g(3) + (m(5) * t16 - t145 * (-pkin(10) * t20 - t16) + t140 * t21 + t142 * t20) * g(2) + (m(5) * t18 - t145 * (-pkin(10) * t24 - t18) + t140 * t25 + t142 * t24) * g(1), t149 * (-g(1) * t24 - g(2) * t20 - g(3) * t46) (-t148 * (t46 * t73 - t62 * t77) + t104 * t11) * g(3) + (-t104 * t146 - t147 * t148) * g(2) + (t104 * t5 - t148 * t6) * g(1) (-g(1) * t5 + g(2) * t146 - g(3) * t11) * m(7)];
taug  = t1(:);
