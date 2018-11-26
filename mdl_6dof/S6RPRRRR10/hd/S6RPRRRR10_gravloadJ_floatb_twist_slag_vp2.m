% Calculate Gravitation load on the joints for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:38:51
% EndTime: 2018-11-23 16:38:53
% DurationCPUTime: 1.58s
% Computational Cost: add. (3613->166), mult. (3737->218), div. (0->0), fcn. (3718->24), ass. (0->93)
t153 = pkin(6) - pkin(13);
t130 = cos(t153) / 0.2e1;
t152 = pkin(6) + pkin(13);
t142 = cos(t152);
t116 = t130 + t142 / 0.2e1;
t157 = sin(pkin(13));
t172 = sin(qJ(1));
t99 = cos(qJ(1));
t113 = -t99 * t116 + t172 * t157;
t90 = sin(pkin(6));
t163 = t90 * t99;
t89 = sin(pkin(7));
t92 = cos(pkin(7));
t56 = t113 * t89 - t92 * t163;
t94 = sin(qJ(4));
t166 = t56 * t94;
t155 = pkin(7) + qJ(3);
t134 = sin(t155) / 0.2e1;
t156 = pkin(7) - qJ(3);
t144 = sin(t156);
t176 = t134 - t144 / 0.2e1;
t129 = sin(t152) / 0.2e1;
t141 = sin(t153);
t74 = t129 - t141 / 0.2e1;
t91 = cos(pkin(13));
t68 = t172 * t91 + t99 * t74;
t135 = cos(t155) / 0.2e1;
t145 = cos(t156);
t76 = t135 - t145 / 0.2e1;
t98 = cos(qJ(3));
t40 = t113 * t176 - t76 * t163 - t68 * t98;
t97 = cos(qJ(4));
t191 = -t40 * t97 + t166;
t185 = mrSges(6,2) - mrSges(7,3);
t93 = sin(qJ(6));
t96 = cos(qJ(6));
t180 = mrSges(7,1) * t96 - mrSges(7,2) * t93 + mrSges(6,1);
t182 = t40 * t94 + t56 * t97;
t88 = qJ(4) + qJ(5);
t85 = sin(t88);
t86 = cos(t88);
t14 = -t40 * t86 + t56 * t85;
t13 = t40 * t85 + t56 * t86;
t186 = m(7) * pkin(5) + t180;
t143 = -m(7) * pkin(12) + t185;
t109 = t172 * t116 + t99 * t157;
t149 = t90 * t172;
t69 = -t172 * t74 + t99 * t91;
t179 = -t109 * t176 - t76 * t149 + t69 * t98;
t184 = t109 * t89 + t92 * t149;
t19 = -t179 * t94 + t184 * t97;
t115 = t129 + t141 / 0.2e1;
t158 = cos(pkin(6));
t75 = t130 - t142 / 0.2e1;
t178 = t115 * t176 - t158 * t76 + t75 * t98;
t67 = -t115 * t89 + t158 * t92;
t183 = -t178 * t94 + t67 * t97;
t177 = -m(6) - m(7);
t131 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t174 = -mrSges(7,1) * t93 - mrSges(7,2) * t96 + t131;
t117 = t134 + t144 / 0.2e1;
t114 = t90 * t117;
t119 = t145 / 0.2e1 + t135;
t95 = sin(qJ(3));
t36 = t113 * t119 + t99 * t114 + t68 * t95;
t173 = m(5) * pkin(3) + t97 * mrSges(5,1) - t94 * mrSges(5,2) - t143 * t85 + t186 * t86 + mrSges(4,1);
t159 = t99 * pkin(1) + qJ(2) * t149;
t148 = t13 * pkin(5) + pkin(12) * t14;
t17 = t179 * t85 - t184 * t86;
t18 = t179 * t86 + t184 * t85;
t147 = -t17 * pkin(5) + pkin(12) * t18;
t26 = -t178 * t85 + t67 * t86;
t27 = t178 * t86 + t67 * t85;
t146 = t26 * pkin(5) + pkin(12) * t27;
t139 = t182 * pkin(4);
t138 = t19 * pkin(4);
t137 = t183 * pkin(4);
t136 = -t172 * pkin(1) + qJ(2) * t163;
t126 = -t180 * t13 + t185 * t14;
t124 = t180 * t17 + t185 * t18;
t121 = -t180 * t26 + t185 * t27;
t107 = -t68 * pkin(2) - t56 * pkin(9) + t136;
t106 = t69 * pkin(2) + t184 * pkin(9) + t159;
t103 = t184 * t94;
t100 = -pkin(11) - pkin(10);
t41 = t109 * t119 - t172 * t114 + t69 * t95;
t84 = pkin(4) * t97 + pkin(3);
t101 = pkin(4) * t103 - t41 * t100 + t179 * t84 + t106;
t46 = -t115 * t119 - t158 * t117 + t75 * t95;
t20 = t179 * t97 + t103;
t2 = t18 * t96 + t41 * t93;
t1 = -t18 * t93 + t41 * t96;
t3 = [(-t99 * mrSges(2,1) + t172 * mrSges(2,2) - m(3) * t159 - t69 * mrSges(3,1) + t109 * mrSges(3,2) - mrSges(3,3) * t149 - m(4) * t106 - t179 * mrSges(4,1) - t184 * mrSges(4,3) - m(5) * (pkin(3) * t179 + t106) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t101 - t18 * mrSges(6,1) - m(7) * (t18 * pkin(5) + t101) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t143 * t17 + t131 * t41) * g(2) + (t172 * mrSges(2,1) + t99 * mrSges(2,2) - m(3) * t136 + t68 * mrSges(3,1) - t113 * mrSges(3,2) - mrSges(3,3) * t163 - m(4) * t107 - t40 * mrSges(4,1) + t56 * mrSges(4,3) - m(5) * (t40 * pkin(3) + t107) + t191 * mrSges(5,1) + t182 * mrSges(5,2) + t143 * t13 + t186 * t14 - t174 * t36 + t177 * (-pkin(4) * t166 + t100 * t36 + t40 * t84 + t107)) * g(1) (-t158 * g(3) + (-t172 * g(1) + t99 * g(2)) * t90) * (m(3) + m(4) + m(5) - t177) (t177 * (-t100 * t178 - t46 * t84) + t174 * t178 + t173 * t46) * g(3) + (t177 * (t100 * t40 - t36 * t84) - t174 * t40 + t173 * t36) * g(2) + (t177 * (-t100 * t179 - t41 * t84) + t174 * t179 + t173 * t41) * g(1) (-t183 * mrSges(5,1) - (-t178 * t97 - t67 * t94) * mrSges(5,2) - m(6) * t137 - m(7) * (t137 + t146) + t121) * g(3) + (-t182 * mrSges(5,1) + t191 * mrSges(5,2) - m(6) * t139 - m(7) * (t139 + t148) + t126) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t138 - m(7) * (t138 + t147) + t124) * g(1) (-m(7) * t146 + t121) * g(3) + (-m(7) * t148 + t126) * g(2) + (-m(7) * t147 + t124) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t93 + t36 * t96) * mrSges(7,1) + (-t14 * t96 - t36 * t93) * mrSges(7,2)) - g(3) * ((-t27 * t93 + t46 * t96) * mrSges(7,1) + (-t27 * t96 - t46 * t93) * mrSges(7,2))];
taug  = t3(:);
