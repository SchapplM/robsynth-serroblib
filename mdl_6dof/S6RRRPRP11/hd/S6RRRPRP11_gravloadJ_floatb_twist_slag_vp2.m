% Calculate Gravitation load on the joints for
% S6RRRPRP11
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
% Datum: 2018-11-23 17:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:48:53
% EndTime: 2018-11-23 17:48:54
% DurationCPUTime: 1.25s
% Computational Cost: add. (1505->138), mult. (1867->176), div. (0->0), fcn. (1860->14), ass. (0->79)
t138 = -m(5) - m(7);
t132 = m(6) - t138;
t137 = mrSges(4,2) - mrSges(5,3);
t145 = m(7) * pkin(5);
t60 = sin(qJ(5));
t71 = -qJ(4) * t132 - t145 * t60 + t137;
t144 = -m(7) * (-qJ(6) - pkin(10)) + mrSges(7,3) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t107 = -mrSges(6,1) - mrSges(7,1);
t136 = mrSges(6,2) + mrSges(7,2);
t61 = sin(qJ(3));
t112 = t60 * t61;
t65 = cos(qJ(3));
t127 = t112 * t145 - t137 * t61 + t144 * t65 + mrSges(3,1);
t123 = m(6) * pkin(10) + t144;
t141 = pkin(3) * t132 + t123;
t63 = sin(qJ(1));
t66 = cos(qJ(2));
t110 = t63 * t66;
t119 = cos(qJ(1));
t100 = pkin(6) + qJ(2);
t88 = sin(t100);
t84 = t88 / 0.2e1;
t101 = pkin(6) - qJ(2);
t89 = sin(t101);
t76 = t84 - t89 / 0.2e1;
t37 = t119 * t76 + t110;
t58 = sin(pkin(6));
t97 = t58 * t119;
t17 = t37 * t61 + t65 * t97;
t62 = sin(qJ(2));
t86 = cos(t100) / 0.2e1;
t90 = cos(t101);
t70 = t90 / 0.2e1 + t86;
t36 = -t119 * t70 + t62 * t63;
t64 = cos(qJ(5));
t140 = -t17 * t60 - t36 * t64;
t139 = t17 * t64 - t36 * t60;
t104 = qJ(4) * t61;
t116 = t36 * t65;
t135 = -pkin(3) * t116 - t36 * t104;
t39 = t119 * t62 + t63 * t70;
t115 = t39 * t65;
t134 = -pkin(3) * t115 - t39 * t104;
t85 = t89 / 0.2e1;
t47 = t84 + t85;
t114 = t47 * t65;
t133 = pkin(3) * t114 + t47 * t104;
t131 = t107 - t145;
t126 = -m(7) * (pkin(5) * t64 + pkin(4)) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t130 = m(6) * pkin(4) + pkin(9) * (m(4) + t132) - t126;
t125 = t107 * t60 - t136 * t64 + t71;
t124 = m(6) * (pkin(4) + pkin(9)) - t126;
t113 = t58 * t63;
t111 = t61 * t64;
t105 = t119 * pkin(1) + pkin(8) * t113;
t103 = cos(pkin(6));
t96 = t119 * t66;
t40 = -t63 * t76 + t96;
t99 = t40 * pkin(2) + t105;
t95 = -pkin(1) * t63 + pkin(8) * t97;
t30 = t36 * pkin(2);
t69 = t85 - t88 / 0.2e1;
t38 = -t119 * t69 + t110;
t94 = pkin(9) * t38 - t30;
t32 = t39 * pkin(2);
t41 = t63 * t69 + t96;
t93 = pkin(9) * t41 - t32;
t46 = t47 * pkin(2);
t48 = t86 - t90 / 0.2e1;
t92 = -pkin(9) * t48 + t46;
t18 = t37 * t65 - t61 * t97;
t87 = -t37 * pkin(2) + t95;
t21 = -t113 * t65 + t40 * t61;
t5 = t21 * t64 - t39 * t60;
t35 = t103 * t61 - t48 * t65;
t34 = -t103 * t65 - t48 * t61;
t22 = t113 * t61 + t40 * t65;
t6 = t21 * t60 + t39 * t64;
t1 = [(-t119 * mrSges(2,1) - m(3) * t105 - t40 * mrSges(3,1) - m(4) * t99 + (-mrSges(3,3) * t58 + mrSges(2,2)) * t63 + t107 * t6 - t136 * t5 + t71 * t21 - t123 * t22 - t130 * t39 - t132 * (t22 * pkin(3) + t99)) * g(2) + (t63 * mrSges(2,1) + t119 * mrSges(2,2) - m(3) * t95 + t37 * mrSges(3,1) - mrSges(3,3) * t97 - m(4) * t87 + t107 * t140 + t136 * t139 - t71 * t17 + t123 * t18 + t130 * t36 + t132 * (pkin(3) * t18 - t87)) * g(1) (-m(4) * t92 - m(6) * (pkin(10) * t114 + t133 + t46) + t138 * (t92 + t133) + t107 * (t112 * t47 - t48 * t64) - t136 * (t111 * t47 + t48 * t60) + t124 * t48 - t127 * t47) * g(3) + (-m(4) * t94 - m(6) * (-pkin(10) * t116 + t135 - t30) + t138 * (t94 + t135) + t107 * (-t112 * t36 + t38 * t64) - t136 * (-t111 * t36 - t38 * t60) - t124 * t38 + t127 * t36) * g(2) + (-m(4) * t93 - m(6) * (-pkin(10) * t115 + t134 - t32) - t136 * (-t111 * t39 - t41 * t60) + t138 * (t93 + t134) + t107 * (-t112 * t39 + t41 * t64) - t124 * t41 + t127 * t39) * g(1) (t125 * t35 + t141 * t34) * g(3) + (t125 * t18 + t141 * t17) * g(2) + (t125 * t22 + t141 * t21) * g(1), t132 * (-g(1) * t21 - g(2) * t17 - g(3) * t34) (-t136 * (-t34 * t60 + t47 * t64) + t131 * (t34 * t64 + t47 * t60)) * g(3) + (t131 * t139 - t136 * t140) * g(2) + (t131 * t5 + t136 * t6) * g(1) (-g(1) * t22 - g(2) * t18 - g(3) * t35) * m(7)];
taug  = t1(:);
