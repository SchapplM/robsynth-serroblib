% Calculate Gravitation load on the joints for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:15:53
% EndTime: 2018-11-23 15:15:55
% DurationCPUTime: 1.44s
% Computational Cost: add. (3307->170), mult. (2889->237), div. (0->0), fcn. (2731->28), ass. (0->96)
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t165 = -m(6) - m(7);
t173 = -t103 * mrSges(7,1) - t107 * mrSges(7,2) + pkin(10) * t165 + mrSges(5,2) - mrSges(6,3);
t172 = m(7) * pkin(5) + t107 * mrSges(7,1) - t103 * mrSges(7,2) + mrSges(6,1);
t139 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t171 = m(5) - t165;
t96 = pkin(7) + qJ(3);
t87 = cos(t96) / 0.2e1;
t97 = pkin(7) - qJ(3);
t94 = cos(t97);
t77 = t94 / 0.2e1 + t87;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t168 = -t139 * t104 + t172 * t108 + mrSges(5,1);
t164 = sin(t97);
t101 = cos(pkin(12));
t106 = sin(qJ(2));
t150 = pkin(6) + qJ(2);
t138 = cos(t150) / 0.2e1;
t151 = pkin(6) - qJ(2);
t143 = cos(t151);
t79 = t143 / 0.2e1 + t138;
t98 = sin(pkin(12));
t54 = t101 * t79 - t106 * t98;
t110 = cos(qJ(2));
t137 = sin(t150) / 0.2e1;
t142 = sin(t151);
t75 = t137 - t142 / 0.2e1;
t55 = t101 * t75 + t110 * t98;
t145 = pkin(3) * t164;
t149 = sin(t96) / 0.2e1;
t82 = pkin(3) * t149;
t99 = sin(pkin(7));
t63 = -t145 / 0.2e1 + t82 - t99 * (pkin(9) + qJ(4));
t109 = cos(qJ(3));
t90 = pkin(3) * t109 + pkin(2);
t163 = t54 * t90 - t55 * t63;
t57 = -t101 * t106 - t79 * t98;
t58 = t101 * t110 - t98 * t75;
t162 = t57 * t90 - t58 * t63;
t74 = t137 + t142 / 0.2e1;
t78 = t138 - t143 / 0.2e1;
t161 = t78 * t63 + t74 * t90;
t100 = sin(pkin(6));
t160 = t100 * t98;
t159 = t104 * t99;
t105 = sin(qJ(3));
t158 = t105 * t55;
t157 = t105 * t58;
t156 = t108 * t99;
t155 = t78 * t105;
t154 = cos(pkin(6));
t153 = t100 * t101;
t102 = cos(pkin(7));
t152 = t100 * t102;
t95 = qJ(3) + pkin(13);
t141 = pkin(7) - t95;
t140 = pkin(7) + t95;
t136 = cos(t140);
t135 = sin(t141);
t70 = t145 / 0.2e1 + t82;
t71 = t77 * pkin(3);
t132 = -pkin(3) * t157 + t70 * t160 + t57 * t71;
t131 = pkin(3) * t155 + t154 * t70 + t74 * t71;
t128 = cos(t141) / 0.2e1;
t127 = sin(t140) / 0.2e1;
t123 = -pkin(3) * t158 - t153 * t70 + t54 * t71;
t122 = -m(4) * pkin(2) - t109 * mrSges(4,1) + t105 * mrSges(4,2) - mrSges(3,1);
t68 = t127 - t135 / 0.2e1;
t69 = t128 - t136 / 0.2e1;
t92 = cos(t95);
t120 = t160 * t69 + t57 * t68 + t58 * t92;
t119 = -t153 * t69 + t54 * t68 + t55 * t92;
t117 = t154 * t69 + t74 * t68 - t78 * t92;
t73 = t149 - t164 / 0.2e1;
t115 = -t73 * mrSges(4,1) - t77 * mrSges(4,2) - mrSges(3,2) + (m(4) * pkin(9) + mrSges(4,3) + mrSges(5,3)) * t99;
t114 = t136 / 0.2e1 + t128;
t113 = t135 / 0.2e1 + t127;
t112 = t100 * t113;
t91 = sin(t95);
t76 = t87 - t94 / 0.2e1;
t72 = t149 + t164 / 0.2e1;
t53 = t102 * t154 - t74 * t99;
t45 = t152 * t98 - t57 * t99;
t44 = -t101 * t152 - t54 * t99;
t33 = t68 * t78 + t74 * t92;
t28 = -t113 * t154 - t114 * t74 - t78 * t91;
t26 = t57 * t92 - t58 * t68;
t24 = t54 * t92 - t55 * t68;
t16 = -t112 * t98 - t114 * t57 + t58 * t91;
t13 = t101 * t112 - t114 * t54 + t55 * t91;
t10 = t104 * t53 + t108 * t117;
t4 = t104 * t45 + t108 * t120;
t2 = t104 * t44 + t108 * t119;
t1 = [(-m(2) - m(3) - m(4) - t171) * g(3) (-m(5) * t161 - t33 * mrSges(5,1) - t172 * (t108 * t33 - t159 * t78) + t173 * (-t114 * t78 + t74 * t91) + t139 * (t104 * t33 + t156 * t78) + t122 * t74 + t115 * t78 + t165 * (t33 * pkin(4) + t161)) * g(3) + (-m(5) * t163 - t24 * mrSges(5,1) + t139 * (t104 * t24 - t156 * t55) - t172 * (t108 * t24 + t159 * t55) + t173 * (t114 * t55 + t54 * t91) + t122 * t54 - t115 * t55 + t165 * (t24 * pkin(4) + t163)) * g(2) + (-m(5) * t162 - t26 * mrSges(5,1) + t139 * (t104 * t26 - t156 * t58) - t172 * (t108 * t26 + t159 * t58) + t173 * (t114 * t58 + t57 * t91) + t122 * t57 - t115 * t58 + t165 * (t26 * pkin(4) + t162)) * g(1) (-(t154 * t72 + t74 * t77 + t155) * mrSges(4,1) - (t78 * t109 + t154 * t76 - t74 * t73) * mrSges(4,2) - m(5) * t131 + t165 * (-t28 * pkin(4) + t131) + t173 * t117 + t168 * t28) * g(3) + (-(-t153 * t72 + t54 * t77 - t158) * mrSges(4,1) - (-t109 * t55 - t153 * t76 - t54 * t73) * mrSges(4,2) - m(5) * t123 + t165 * (-t13 * pkin(4) + t123) + t173 * t119 + t168 * t13) * g(2) + (-(t72 * t160 + t57 * t77 - t157) * mrSges(4,1) - (-t109 * t58 + t76 * t160 - t57 * t73) * mrSges(4,2) - m(5) * t132 + t165 * (-t16 * pkin(4) + t132) + t173 * t120 + t168 * t16) * g(1), t171 * (-g(1) * t45 - g(2) * t44 - g(3) * t53) (-t172 * (-t104 * t117 + t108 * t53) + t139 * t10) * g(3) + (t139 * t2 - t172 * (-t104 * t119 + t108 * t44)) * g(2) + (t139 * t4 - t172 * (-t104 * t120 + t108 * t45)) * g(1), -g(1) * ((-t103 * t4 + t107 * t16) * mrSges(7,1) + (-t103 * t16 - t107 * t4) * mrSges(7,2)) - g(2) * ((-t103 * t2 + t107 * t13) * mrSges(7,1) + (-t103 * t13 - t107 * t2) * mrSges(7,2)) - g(3) * ((-t10 * t103 + t107 * t28) * mrSges(7,1) + (-t10 * t107 - t103 * t28) * mrSges(7,2))];
taug  = t1(:);
