% Calculate Gravitation load on the joints for
% S6RPRPRR11
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
% Datum: 2018-11-23 16:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:09:22
% EndTime: 2018-11-23 16:09:23
% DurationCPUTime: 1.30s
% Computational Cost: add. (2992->143), mult. (3087->191), div. (0->0), fcn. (3041->24), ass. (0->77)
t71 = sin(pkin(6));
t80 = cos(qJ(1));
t134 = t71 * t80;
t126 = pkin(7) + qJ(3);
t111 = sin(t126) / 0.2e1;
t127 = pkin(7) - qJ(3);
t118 = sin(t127);
t141 = t111 - t118 / 0.2e1;
t136 = sin(qJ(1));
t123 = pkin(6) + pkin(12);
t107 = sin(t123) / 0.2e1;
t124 = pkin(6) - pkin(12);
t115 = sin(t124);
t54 = t107 - t115 / 0.2e1;
t73 = cos(pkin(12));
t48 = t136 * t73 + t80 * t54;
t112 = cos(t126) / 0.2e1;
t119 = cos(t127);
t56 = t112 - t119 / 0.2e1;
t79 = cos(qJ(3));
t128 = sin(pkin(12));
t108 = cos(t124) / 0.2e1;
t116 = cos(t123);
t95 = t108 + t116 / 0.2e1;
t92 = t136 * t128 - t80 * t95;
t23 = -t56 * t134 + t141 * t92 - t48 * t79;
t70 = sin(pkin(7));
t74 = cos(pkin(7));
t37 = -t74 * t134 + t92 * t70;
t68 = pkin(13) + qJ(5);
t65 = sin(t68);
t66 = cos(t68);
t4 = -t23 * t66 + t37 * t65;
t148 = t23 * t65 + t37 * t66;
t76 = sin(qJ(6));
t78 = cos(qJ(6));
t142 = m(7) * pkin(5) + mrSges(7,1) * t78 - mrSges(7,2) * t76 + mrSges(6,1);
t117 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t121 = t71 * t136;
t88 = t80 * t128 + t136 * t95;
t84 = -t74 * t121 - t88 * t70;
t49 = -t136 * t54 + t80 * t73;
t145 = -t56 * t121 - t141 * t88 + t49 * t79;
t129 = cos(pkin(6));
t55 = t108 - t116 / 0.2e1;
t94 = t107 + t115 / 0.2e1;
t144 = -t129 * t56 + t141 * t94 + t55 * t79;
t143 = -m(6) - m(7);
t106 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t137 = -mrSges(7,1) * t76 - mrSges(7,2) * t78 + t106;
t140 = m(5) - t143;
t77 = sin(qJ(3));
t96 = t111 + t118 / 0.2e1;
t93 = t71 * t96;
t98 = t119 / 0.2e1 + t112;
t19 = t48 * t77 + t80 * t93 + t92 * t98;
t69 = sin(pkin(13));
t72 = cos(pkin(13));
t138 = m(5) * pkin(3) + t72 * mrSges(5,1) - t69 * mrSges(5,2) - t117 * t65 + t142 * t66 + mrSges(4,1);
t135 = t37 * t69;
t130 = t80 * pkin(1) + qJ(2) * t121;
t113 = -t136 * pkin(1) + qJ(2) * t134;
t86 = -t48 * pkin(2) - t37 * pkin(9) + t113;
t85 = t49 * pkin(2) - t84 * pkin(9) + t130;
t82 = t84 * t69;
t24 = -t136 * t93 + t49 * t77 + t88 * t98;
t64 = pkin(4) * t72 + pkin(3);
t75 = -pkin(10) - qJ(4);
t81 = -pkin(4) * t82 + t145 * t64 - t24 * t75 + t85;
t47 = t129 * t74 - t94 * t70;
t29 = -t129 * t96 + t55 * t77 - t94 * t98;
t10 = t144 * t66 + t47 * t65;
t8 = t145 * t66 - t84 * t65;
t7 = t145 * t65 + t84 * t66;
t2 = t24 * t76 + t78 * t8;
t1 = t24 * t78 - t76 * t8;
t3 = [(-t80 * mrSges(2,1) + t136 * mrSges(2,2) - m(3) * t130 - t49 * mrSges(3,1) + t88 * mrSges(3,2) - mrSges(3,3) * t121 - m(4) * t85 - t145 * mrSges(4,1) + t84 * mrSges(4,3) - m(5) * (pkin(3) * t145 + t85) - (t145 * t72 - t82) * mrSges(5,1) - (-t145 * t69 - t84 * t72) * mrSges(5,2) - m(6) * t81 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t81) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t117 * t7 + t106 * t24) * g(2) + (t136 * mrSges(2,1) + t80 * mrSges(2,2) - m(3) * t113 + t48 * mrSges(3,1) - t92 * mrSges(3,2) - mrSges(3,3) * t134 - m(4) * t86 - t23 * mrSges(4,1) + t37 * mrSges(4,3) - m(5) * (t23 * pkin(3) + t86) - (t23 * t72 - t135) * mrSges(5,1) - (-t23 * t69 - t37 * t72) * mrSges(5,2) + t117 * t148 + t142 * t4 - t137 * t19 + t143 * (-pkin(4) * t135 + t19 * t75 + t23 * t64 + t86)) * g(1) (-t129 * g(3) + (-g(1) * t136 + g(2) * t80) * t71) * (m(3) + m(4) + t140) (t143 * (-t144 * t75 - t29 * t64) + t137 * t144 + t138 * t29) * g(3) + (t143 * (-t19 * t64 + t23 * t75) - t137 * t23 + t138 * t19) * g(2) + (t143 * (-t145 * t75 - t24 * t64) + t137 * t145 + t138 * t24) * g(1), t140 * (-g(1) * t24 - g(2) * t19 - g(3) * t29) (-t142 * (-t144 * t65 + t47 * t66) + t117 * t10) * g(3) + (t117 * t4 - t142 * t148) * g(2) + (t117 * t8 + t142 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t78 - t4 * t76) * mrSges(7,1) + (-t19 * t76 - t4 * t78) * mrSges(7,2)) - g(3) * ((-t10 * t76 + t29 * t78) * mrSges(7,1) + (-t10 * t78 - t29 * t76) * mrSges(7,2))];
taug  = t3(:);
