% Calculate Gravitation load on the joints for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:34:29
% EndTime: 2018-11-23 15:34:31
% DurationCPUTime: 1.34s
% Computational Cost: add. (3715->156), mult. (3716->222), div. (0->0), fcn. (3624->24), ass. (0->93)
t179 = mrSges(6,2) - mrSges(7,3);
t95 = sin(qJ(6));
t98 = cos(qJ(6));
t174 = -t98 * mrSges(7,1) + t95 * mrSges(7,2) - mrSges(6,1);
t180 = m(7) * pkin(5) - t174;
t135 = -m(7) * pkin(12) + t179;
t100 = cos(qJ(3));
t148 = pkin(6) + qJ(2);
t128 = sin(t148) / 0.2e1;
t149 = pkin(6) - qJ(2);
t137 = sin(t149);
t113 = t128 + t137 / 0.2e1;
t153 = cos(pkin(6));
t146 = pkin(7) + qJ(3);
t127 = sin(t146) / 0.2e1;
t147 = pkin(7) - qJ(3);
t136 = sin(t147);
t81 = t127 - t136 / 0.2e1;
t129 = cos(t146) / 0.2e1;
t138 = cos(t147);
t83 = t129 - t138 / 0.2e1;
t130 = cos(t148) / 0.2e1;
t139 = cos(t149);
t84 = t130 - t139 / 0.2e1;
t107 = -t84 * t100 + t113 * t81 - t153 * t83;
t93 = sin(pkin(7));
t94 = cos(pkin(7));
t71 = -t113 * t93 + t153 * t94;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t177 = -t107 * t96 + t71 * t99;
t115 = t139 / 0.2e1 + t130;
t150 = sin(pkin(13));
t152 = cos(pkin(13));
t170 = sin(qJ(2));
t106 = t150 * t115 + t152 * t170;
t151 = sin(pkin(6));
t123 = t151 * t150;
t101 = cos(qJ(2));
t82 = t128 - t137 / 0.2e1;
t74 = t152 * t101 - t150 * t82;
t103 = t74 * t100 - t106 * t81 - t83 * t123;
t58 = t106 * t93 + t94 * t123;
t176 = -t103 * t96 + t58 * t99;
t105 = -t152 * t115 + t150 * t170;
t124 = t152 * t151;
t72 = t150 * t101 + t152 * t82;
t104 = t72 * t100 - t105 * t81 + t83 * t124;
t57 = t105 * t93 - t94 * t124;
t175 = -t104 * t96 + t57 * t99;
t171 = m(6) + m(7);
t145 = m(4) + m(5) + t171;
t173 = pkin(2) * t145 + mrSges(3,1);
t117 = -m(5) * pkin(3) - t99 * mrSges(5,1) + t96 * mrSges(5,2) - mrSges(4,1);
t111 = -m(5) * pkin(10) - t95 * mrSges(7,1) - t98 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t92 = qJ(4) + qJ(5);
t90 = sin(t92);
t91 = cos(t92);
t172 = -t135 * t90 + t180 * t91 - t117;
t163 = t90 * t93;
t162 = t91 * t93;
t11 = -t104 * t90 + t57 * t91;
t12 = t104 * t91 + t57 * t90;
t142 = t11 * pkin(5) + pkin(12) * t12;
t13 = -t103 * t90 + t58 * t91;
t14 = t103 * t91 + t58 * t90;
t141 = t13 * pkin(5) + pkin(12) * t14;
t24 = -t107 * t90 + t71 * t91;
t25 = t107 * t91 + t71 * t90;
t140 = t24 * pkin(5) + pkin(12) * t25;
t134 = t175 * pkin(4);
t133 = t176 * pkin(4);
t132 = t177 * pkin(4);
t122 = t174 * t11 + t179 * t12;
t121 = t174 * t13 + t179 * t14;
t118 = t174 * t24 + t179 * t25;
t114 = t138 / 0.2e1 + t129;
t112 = t127 + t136 / 0.2e1;
t110 = t112 * t151;
t108 = -mrSges(3,2) + (mrSges(5,2) * t99 + mrSges(4,3) + (t171 * pkin(4) + mrSges(5,1)) * t96 + t145 * pkin(9)) * t93;
t102 = -pkin(11) - pkin(10);
t97 = sin(qJ(3));
t89 = pkin(4) * t99 + pkin(3);
t56 = t113 * t100 + t84 * t81;
t55 = t113 * t97 - t84 * t114;
t48 = -t153 * t112 - t113 * t114 - t84 * t97;
t46 = -t106 * t100 - t74 * t81;
t45 = -t106 * t97 + t114 * t74;
t44 = -t105 * t100 - t72 * t81;
t43 = -t105 * t97 + t114 * t72;
t33 = t106 * t114 - t150 * t110 + t74 * t97;
t30 = t105 * t114 + t152 * t110 + t72 * t97;
t1 = [(-m(2) - m(3) - t145) * g(3) (t135 * (t84 * t162 + t56 * t90) + t117 * t56 - t180 * (-t84 * t163 + t56 * t91) + t111 * t55 + t108 * t84 - t171 * (-t55 * t102 + t56 * t89) - t173 * t113) * g(3) + (t135 * (-t162 * t72 + t44 * t90) + t117 * t44 - t180 * (t163 * t72 + t44 * t91) + t111 * t43 - t108 * t72 - t171 * (-t43 * t102 + t44 * t89) + t173 * t105) * g(2) + (t135 * (-t162 * t74 + t46 * t90) + t117 * t46 - t180 * (t163 * t74 + t46 * t91) + t111 * t45 - t108 * t74 - t171 * (-t45 * t102 + t46 * t89) + t173 * t106) * g(1) (-t171 * (-t102 * t107 - t48 * t89) + t111 * t107 + t172 * t48) * g(3) + (-t171 * (-t102 * t104 - t30 * t89) + t111 * t104 + t172 * t30) * g(2) + (-t171 * (-t102 * t103 - t33 * t89) + t111 * t103 + t172 * t33) * g(1) (-t177 * mrSges(5,1) - (-t107 * t99 - t71 * t96) * mrSges(5,2) - m(6) * t132 - m(7) * (t132 + t140) + t118) * g(3) + (-t175 * mrSges(5,1) - (-t104 * t99 - t57 * t96) * mrSges(5,2) - m(6) * t134 - m(7) * (t134 + t142) + t122) * g(2) + (-t176 * mrSges(5,1) - (-t103 * t99 - t58 * t96) * mrSges(5,2) - m(6) * t133 - m(7) * (t133 + t141) + t121) * g(1) (-m(7) * t140 + t118) * g(3) + (-m(7) * t142 + t122) * g(2) + (-m(7) * t141 + t121) * g(1), -g(1) * ((-t14 * t95 + t33 * t98) * mrSges(7,1) + (-t14 * t98 - t33 * t95) * mrSges(7,2)) - g(2) * ((-t12 * t95 + t30 * t98) * mrSges(7,1) + (-t12 * t98 - t30 * t95) * mrSges(7,2)) - g(3) * ((-t25 * t95 + t48 * t98) * mrSges(7,1) + (-t25 * t98 - t48 * t95) * mrSges(7,2))];
taug  = t1(:);
