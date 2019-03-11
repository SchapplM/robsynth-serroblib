% Calculate Gravitation load on the joints for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:35
% EndTime: 2019-03-09 20:24:39
% DurationCPUTime: 1.58s
% Computational Cost: add. (1219->153), mult. (3291->228), div. (0->0), fcn. (4125->14), ass. (0->78)
t146 = m(6) + m(7);
t142 = m(5) + t146;
t100 = -qJ(4) * t142 + mrSges(4,2) - mrSges(5,3);
t101 = -pkin(11) * t146 - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t78 = sin(qJ(6));
t82 = cos(qJ(6));
t89 = -mrSges(7,1) * t78 - mrSges(7,2) * t82 + t101;
t145 = m(7) * pkin(5) + mrSges(7,1) * t82 - mrSges(7,2) * t78 + mrSges(6,1);
t150 = m(7) * pkin(12) - mrSges(6,2) + mrSges(7,3);
t121 = cos(pkin(7));
t77 = sin(pkin(6));
t113 = t77 * t121;
t135 = sin(qJ(1));
t122 = cos(pkin(6));
t105 = t122 * t135;
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t85 = cos(qJ(1));
t59 = -t105 * t84 - t81 * t85;
t76 = sin(pkin(7));
t42 = t113 * t135 - t59 * t76;
t111 = t85 * t122;
t57 = -t111 * t84 + t135 * t81;
t149 = t113 * t85 - t57 * t76;
t148 = pkin(3) * t142 - t89;
t136 = cos(qJ(3));
t104 = t121 * t136;
t116 = t77 * t136;
t58 = t111 * t81 + t135 * t84;
t80 = sin(qJ(3));
t21 = t116 * t76 * t85 + t104 * t57 + t58 * t80;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t147 = t149 * t79 + t21 * t83;
t6 = -t149 * t83 + t21 * t79;
t138 = -t145 * t79 + t150 * t83 + t100;
t132 = t76 * t79;
t131 = t76 * t83;
t130 = t77 * t81;
t129 = t77 * t84;
t128 = t77 * t85;
t127 = -mrSges(4,3) - mrSges(5,1);
t119 = t76 * t130;
t124 = pkin(2) * t129 + pkin(10) * t119;
t115 = t77 * t135;
t123 = pkin(1) * t85 + pkin(9) * t115;
t112 = t80 * t121;
t51 = -t112 * t130 + t116 * t84;
t118 = pkin(3) * t51 + t124;
t114 = t76 * t122;
t108 = t76 * t115;
t106 = -pkin(1) * t135 + pkin(9) * t128;
t60 = -t105 * t81 + t84 * t85;
t98 = t60 * pkin(2) + pkin(10) * t42 + t123;
t26 = t60 * t136 + (t121 * t59 + t108) * t80;
t94 = pkin(3) * t26 + t98;
t92 = pkin(4) * t42 + t94;
t91 = -t58 * pkin(2) + pkin(10) * t149 + t106;
t22 = -t128 * t76 * t80 - t112 * t57 + t136 * t58;
t90 = -pkin(3) * t22 + t91;
t87 = mrSges(3,2) + (-t146 * pkin(4) + (-m(4) - t142) * pkin(10) + t127) * t76;
t56 = t121 * t122 - t129 * t76;
t54 = t59 * pkin(2);
t52 = t57 * pkin(2);
t50 = (t104 * t81 + t80 * t84) * t77;
t37 = t80 * t114 + (t112 * t84 + t136 * t81) * t77;
t36 = -t104 * t129 - t114 * t136 + t130 * t80;
t32 = -t112 * t60 + t136 * t59;
t31 = t104 * t60 + t59 * t80;
t30 = -t112 * t58 - t136 * t57;
t29 = t104 * t58 - t57 * t80;
t25 = -t104 * t59 - t108 * t136 + t60 * t80;
t20 = t36 * t79 + t56 * t83;
t8 = t25 * t79 + t42 * t83;
t7 = -t25 * t83 + t42 * t79;
t2 = t26 * t78 + t8 * t82;
t1 = t26 * t82 - t78 * t8;
t3 = [(-t85 * mrSges(2,1) + t135 * mrSges(2,2) - m(3) * t123 - t60 * mrSges(3,1) - t59 * mrSges(3,2) - mrSges(3,3) * t115 - m(4) * t98 - m(5) * t94 - m(6) * t92 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t92) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t150 * t7 + t127 * t42 + t100 * t25 + t101 * t26) * g(2) + (t135 * mrSges(2,1) - m(3) * t106 + t58 * mrSges(3,1) - t57 * mrSges(3,2) - m(4) * t91 - m(5) * t90 + (-mrSges(3,3) * t77 + mrSges(2,2)) * t85 + t127 * t149 - t150 * t147 - t100 * t21 + t145 * t6 - t89 * t22 + t146 * (-pkin(4) * t149 - t90)) * g(1) (-m(4) * t124 - m(5) * t118 + t100 * t50 - t150 * (t119 * t79 - t50 * t83) - t145 * (t119 * t83 + t50 * t79) + t89 * t51 + (-mrSges(3,1) * t84 + (t127 * t76 + mrSges(3,2)) * t81) * t77 - t146 * (pkin(4) * t119 + t118)) * g(3) + (mrSges(3,1) * t57 + m(4) * t52 + t150 * (-t132 * t58 + t29 * t83) + t100 * t29 - t145 * (t131 * t58 + t29 * t79) + t89 * t30 + t87 * t58 - t142 * (pkin(3) * t30 - t52)) * g(2) + (-mrSges(3,1) * t59 - m(4) * t54 + t100 * t31 + t150 * (-t132 * t60 + t31 * t83) - t145 * (t131 * t60 + t31 * t79) + t89 * t32 + t87 * t60 - t142 * (pkin(3) * t32 + t54)) * g(1) (t138 * t37 + t148 * t36) * g(3) + (t138 * t22 + t148 * t21) * g(2) + (t138 * t26 + t148 * t25) * g(1), t142 * (-g(1) * t25 - g(2) * t21 - g(3) * t36) (-t150 * t20 - t145 * (t36 * t83 - t56 * t79)) * g(3) + (-t145 * t147 - t150 * t6) * g(2) + (t145 * t7 - t150 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t22 * t82 - t6 * t78) * mrSges(7,1) + (-t22 * t78 - t6 * t82) * mrSges(7,2)) - g(3) * ((-t20 * t78 + t37 * t82) * mrSges(7,1) + (-t20 * t82 - t37 * t78) * mrSges(7,2))];
taug  = t3(:);
