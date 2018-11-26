% Calculate Gravitation load on the joints for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:53:04
% EndTime: 2018-11-23 14:53:06
% DurationCPUTime: 1.30s
% Computational Cost: add. (6696->141), mult. (6658->208), div. (0->0), fcn. (6532->30), ass. (0->99)
t154 = m(6) + m(7);
t133 = m(5) + t154;
t160 = m(3) + m(4) + t133;
t70 = sin(qJ(6));
t73 = cos(qJ(6));
t159 = -t70 * mrSges(7,1) - t73 * mrSges(7,2) - pkin(11) * t154 + mrSges(5,2) - mrSges(6,3);
t158 = m(7) * pkin(5) + t73 * mrSges(7,1) - t70 * mrSges(7,2) + mrSges(6,1);
t120 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t156 = pkin(4) * t154 - t120 * t71 + t158 * t74 + mrSges(5,1);
t147 = sin(qJ(3));
t68 = sin(pkin(8));
t146 = t68 * t71;
t145 = t68 * t74;
t141 = cos(pkin(6));
t140 = cos(pkin(7));
t139 = cos(pkin(13));
t138 = cos(pkin(14));
t137 = sin(pkin(6));
t136 = sin(pkin(7));
t135 = sin(pkin(13));
t134 = sin(pkin(14));
t132 = pkin(7) - qJ(3);
t131 = pkin(7) + qJ(3);
t130 = pkin(8) - qJ(4);
t129 = pkin(8) + qJ(4);
t128 = pkin(6) - pkin(14);
t127 = pkin(6) + pkin(14);
t125 = cos(t132);
t124 = cos(t130);
t123 = cos(t129);
t122 = sin(t132);
t121 = sin(t130);
t119 = cos(t127);
t118 = sin(t128);
t117 = cos(t131) / 0.2e1;
t116 = t124 / 0.2e1;
t115 = t123 / 0.2e1;
t114 = sin(t131) / 0.2e1;
t113 = sin(t129) / 0.2e1;
t110 = t139 * t137;
t109 = t137 * t135;
t108 = cos(t128) / 0.2e1;
t107 = sin(t127) / 0.2e1;
t106 = t116 - t123 / 0.2e1;
t102 = -mrSges(4,2) + (t133 * pkin(10) + mrSges(5,3)) * t68;
t101 = t125 / 0.2e1 + t117;
t100 = t116 + t115;
t99 = t115 - t124 / 0.2e1;
t98 = t114 + t122 / 0.2e1;
t97 = t113 + t121 / 0.2e1;
t95 = t108 - t119 / 0.2e1;
t94 = t108 + t119 / 0.2e1;
t93 = t107 - t118 / 0.2e1;
t92 = t107 + t118 / 0.2e1;
t91 = t98 * t137;
t89 = t92 * t136 - t141 * t140;
t88 = -t135 * t93 + t139 * t138;
t87 = t139 * t134 + t135 * t94;
t86 = t135 * t138 + t139 * t93;
t85 = t135 * t134 - t139 * t94;
t62 = t114 - t122 / 0.2e1;
t63 = t117 - t125 / 0.2e1;
t76 = cos(qJ(3));
t53 = -t141 * t63 + t92 * t62 + t95 * t76;
t84 = -t140 * t109 - t87 * t136;
t83 = t140 * t110 - t85 * t136;
t82 = -t92 * t101 - t141 * t98 + t95 * t147;
t47 = -t63 * t109 - t87 * t62 + t88 * t76;
t45 = t63 * t110 - t85 * t62 + t86 * t76;
t61 = t113 - t121 / 0.2e1;
t75 = cos(qJ(4));
t81 = t53 * t75 - t82 * t61;
t80 = t87 * t101 - t135 * t91 + t88 * t147;
t79 = t85 * t101 + t139 * t91 + t86 * t147;
t78 = t45 * t75 - t79 * t61;
t77 = t47 * t75 - t80 * t61;
t72 = sin(qJ(4));
t69 = cos(pkin(8));
t52 = t82 * pkin(3);
t44 = t80 * pkin(3);
t43 = t79 * pkin(3);
t42 = t82 * t68 - t89 * t69;
t35 = t80 * t68 - t84 * t69;
t34 = t79 * t68 - t83 * t69;
t33 = -t53 * t61 - t82 * t75;
t29 = -t89 * t106 + t81;
t28 = t82 * t100 + t53 * t72 + t89 * t97;
t26 = -t47 * t61 - t80 * t75;
t24 = -t45 * t61 - t79 * t75;
t17 = -t84 * t106 + t77;
t16 = t80 * t100 + t47 * t72 + t84 * t97;
t14 = -t83 * t106 + t78;
t13 = t79 * t100 + t45 * t72 + t83 * t97;
t10 = t29 * t74 + t42 * t71;
t4 = t17 * t74 + t35 * t71;
t2 = t14 * t74 + t34 * t71;
t1 = [(-m(2) - t160) * g(3) (-t109 * g(1) + t110 * g(2) - t141 * g(3)) * t160 (t82 * mrSges(4,1) + m(5) * t52 - t33 * mrSges(5,1) - t102 * t53 - t158 * (t146 * t53 + t33 * t74) + t159 * (t100 * t53 - t82 * t72) + t120 * (-t145 * t53 + t33 * t71) - t154 * (t33 * pkin(4) - t52)) * g(3) + (t79 * mrSges(4,1) + m(5) * t43 - t24 * mrSges(5,1) + t120 * (-t145 * t45 + t24 * t71) - t102 * t45 - t158 * (t146 * t45 + t24 * t74) + t159 * (t100 * t45 - t79 * t72) - t154 * (t24 * pkin(4) - t43)) * g(2) + (t80 * mrSges(4,1) + m(5) * t44 - t26 * mrSges(5,1) + t120 * (-t145 * t47 + t26 * t71) - t102 * t47 - t158 * (t146 * t47 + t26 * t74) + t159 * (t100 * t47 - t80 * t72) - t154 * (t26 * pkin(4) - t44)) * g(1) (t159 * (t89 * t99 + t81) + t156 * t28) * g(3) + (t159 * (t83 * t99 + t78) + t156 * t13) * g(2) + (t159 * (t84 * t99 + t77) + t156 * t16) * g(1) (-t158 * (-t29 * t71 + t42 * t74) + t120 * t10) * g(3) + (t120 * t2 - t158 * (-t14 * t71 + t34 * t74)) * g(2) + (t120 * t4 - t158 * (-t17 * t71 + t35 * t74)) * g(1), -g(1) * ((t16 * t73 - t4 * t70) * mrSges(7,1) + (-t16 * t70 - t4 * t73) * mrSges(7,2)) - g(2) * ((t13 * t73 - t2 * t70) * mrSges(7,1) + (-t13 * t70 - t2 * t73) * mrSges(7,2)) - g(3) * ((-t10 * t70 + t28 * t73) * mrSges(7,1) + (-t10 * t73 - t28 * t70) * mrSges(7,2))];
taug  = t1(:);
