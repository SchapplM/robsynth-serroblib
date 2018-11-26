% Calculate Gravitation load on the joints for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:08:10
% EndTime: 2018-11-23 16:08:11
% DurationCPUTime: 1.39s
% Computational Cost: add. (3266->176), mult. (2928->230), div. (0->0), fcn. (2824->28), ass. (0->107)
t170 = -m(6) - m(7);
t95 = sin(qJ(6));
t99 = cos(qJ(6));
t168 = m(7) * pkin(5) + t99 * mrSges(7,1) - t95 * mrSges(7,2) + mrSges(6,1);
t163 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t169 = mrSges(5,2) - mrSges(6,3);
t172 = -t95 * mrSges(7,1) - t99 * mrSges(7,2) + t169;
t102 = cos(qJ(1));
t90 = sin(pkin(6));
t151 = t102 * t90;
t86 = qJ(3) + pkin(13);
t141 = pkin(7) + t86;
t122 = sin(t141) / 0.2e1;
t142 = pkin(7) - t86;
t131 = sin(t142);
t167 = t122 - t131 / 0.2e1;
t147 = pkin(6) - pkin(12);
t130 = cos(t147) / 0.2e1;
t146 = pkin(6) + pkin(12);
t139 = cos(t146);
t111 = t130 + t139 / 0.2e1;
t150 = sin(pkin(12));
t98 = sin(qJ(1));
t50 = -t102 * t111 + t150 * t98;
t129 = sin(t146) / 0.2e1;
t138 = sin(t147);
t64 = t129 - t138 / 0.2e1;
t91 = cos(pkin(12));
t51 = t102 * t64 + t98 * t91;
t123 = cos(t142) / 0.2e1;
t132 = cos(t141);
t61 = t123 - t132 / 0.2e1;
t83 = cos(t86);
t16 = -t151 * t61 - t167 * t50 + t51 * t83;
t104 = t102 * t150 + t111 * t98;
t153 = t90 * t98;
t52 = t102 * t91 - t98 * t64;
t21 = -t104 * t167 + t153 * t61 + t52 * t83;
t110 = t129 + t138 / 0.2e1;
t65 = t130 - t139 / 0.2e1;
t93 = cos(pkin(6));
t26 = t110 * t167 + t93 * t61 + t65 * t83;
t100 = cos(qJ(5));
t89 = sin(pkin(7));
t92 = cos(pkin(7));
t39 = t151 * t92 - t50 * t89;
t96 = sin(qJ(5));
t171 = -t39 * t100 - t16 * t96;
t4 = t100 * t16 - t39 * t96;
t166 = m(5) - t170;
t165 = t168 * t100 - t163 * t96 + mrSges(5,1);
t164 = -m(4) * pkin(9) - mrSges(4,3) - mrSges(5,3);
t161 = t170 * pkin(10) + t172;
t87 = pkin(7) + qJ(3);
t84 = cos(t87);
t160 = t84 / 0.2e1;
t88 = pkin(7) - qJ(3);
t159 = cos(t88);
t158 = sin(t88);
t157 = pkin(3) * t84;
t97 = sin(qJ(3));
t156 = t51 * t97;
t155 = t52 * t97;
t154 = t65 * t97;
t152 = t102 * pkin(1) + qJ(2) * t153;
t149 = t159 / 0.2e1;
t148 = sin(t87) / 0.2e1;
t145 = pkin(3) * t158;
t143 = -t98 * pkin(1) + qJ(2) * t151;
t72 = pkin(3) * t148;
t94 = pkin(9) + qJ(4);
t53 = -t145 / 0.2e1 + t72 - t89 * t94;
t73 = pkin(3) * t149;
t54 = t73 - t157 / 0.2e1 + t92 * t94;
t101 = cos(qJ(3));
t81 = pkin(3) * t101 + pkin(2);
t137 = -t104 * t53 + t54 * t153 + t52 * t81 + t152;
t62 = t145 / 0.2e1 + t72;
t63 = t73 + t157 / 0.2e1;
t134 = -pkin(3) * t155 - t104 * t63 + t62 * t153;
t133 = -pkin(3) * t154 + t110 * t63 + t93 * t62;
t124 = t54 * t151 + t50 * t53 - t51 * t81 + t143;
t117 = -pkin(3) * t156 - t151 * t62 - t50 * t63;
t108 = t131 / 0.2e1 + t122;
t106 = t90 * t108;
t109 = t132 / 0.2e1 + t123;
t82 = sin(t86);
t20 = t104 * t109 - t106 * t98 + t52 * t82;
t116 = t21 * pkin(4) + pkin(10) * t20 + t137;
t66 = t148 + t158 / 0.2e1;
t69 = t149 + t160;
t115 = t151 * t66 + t50 * t69 + t156;
t67 = t148 - t158 / 0.2e1;
t68 = t160 - t159 / 0.2e1;
t114 = -t51 * t101 - t151 * t68 + t50 * t67;
t41 = t104 * t89 + t153 * t92;
t103 = t52 * t101 - t104 * t67 - t153 * t68;
t15 = t102 * t106 + t109 * t50 + t51 * t82;
t49 = -t110 * t89 + t93 * t92;
t25 = -t108 * t93 - t109 * t110 + t65 * t82;
t23 = -t104 * t69 + t153 * t66 - t155;
t10 = t100 * t26 + t49 * t96;
t8 = t100 * t21 + t41 * t96;
t7 = -t41 * t100 + t21 * t96;
t2 = t20 * t95 + t8 * t99;
t1 = t20 * t99 - t8 * t95;
t3 = [(-mrSges(2,1) * t102 + t98 * mrSges(2,2) - m(3) * t152 - t52 * mrSges(3,1) + t104 * mrSges(3,2) - mrSges(3,3) * t153 - m(4) * (t52 * pkin(2) + t152) - t103 * mrSges(4,1) - t23 * mrSges(4,2) - m(5) * t137 - t21 * mrSges(5,1) - m(6) * t116 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t116) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t163 * t7 + t164 * t41 + t169 * t20) * g(2) + (t98 * mrSges(2,1) + mrSges(2,2) * t102 - m(3) * t143 + t51 * mrSges(3,1) - t50 * mrSges(3,2) - mrSges(3,3) * t151 - m(4) * (-t51 * pkin(2) + t143) - t114 * mrSges(4,1) - t115 * mrSges(4,2) - m(5) * t124 + t16 * mrSges(5,1) + t170 * (-pkin(4) * t16 - pkin(10) * t15 + t124) + t163 * t171 + t164 * t39 + t168 * t4 - t172 * t15) * g(1) (-t93 * g(3) + (-g(1) * t98 + g(2) * t102) * t90) * (m(3) + m(4) + t166) (-(t110 * t69 + t93 * t66 - t154) * mrSges(4,1) - (-t65 * t101 - t110 * t67 + t93 * t68) * mrSges(4,2) - m(5) * t133 + t170 * (-t25 * pkin(4) + t133) + t161 * t26 + t165 * t25) * g(3) + (-m(5) * t117 + t115 * mrSges(4,1) - t114 * mrSges(4,2) + t170 * (-t15 * pkin(4) + t117) + t161 * t16 + t165 * t15) * g(2) + (-m(5) * t134 - t23 * mrSges(4,1) + t103 * mrSges(4,2) + t170 * (-t20 * pkin(4) + t134) + t161 * t21 + t165 * t20) * g(1), t166 * (-g(1) * t41 + g(2) * t39 - g(3) * t49) (-t168 * (t100 * t49 - t26 * t96) + t163 * t10) * g(3) + (t163 * t4 - t168 * t171) * g(2) + (t163 * t8 + t168 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t15 * t99 - t4 * t95) * mrSges(7,1) + (-t15 * t95 - t4 * t99) * mrSges(7,2)) - g(3) * ((-t10 * t95 + t25 * t99) * mrSges(7,1) + (-t10 * t99 - t25 * t95) * mrSges(7,2))];
taug  = t3(:);
