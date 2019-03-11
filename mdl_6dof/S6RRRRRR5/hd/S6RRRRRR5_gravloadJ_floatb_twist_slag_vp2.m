% Calculate Gravitation load on the joints for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:40
% EndTime: 2019-03-10 03:56:44
% DurationCPUTime: 1.47s
% Computational Cost: add. (1080->166), mult. (1442->216), div. (0->0), fcn. (1646->14), ass. (0->81)
t163 = mrSges(6,2) - mrSges(7,3);
t82 = sin(qJ(6));
t86 = cos(qJ(6));
t162 = mrSges(7,1) * t86 - mrSges(7,2) * t82 + mrSges(6,1);
t164 = m(7) * pkin(5) + t162;
t110 = -m(7) * pkin(12) + t163;
t90 = -pkin(10) - pkin(9);
t92 = -m(4) * pkin(9) + m(5) * t90 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t157 = mrSges(7,1) * t82 + mrSges(7,2) * t86 - t92;
t159 = m(6) + m(7);
t124 = cos(pkin(6));
t81 = sin(pkin(6));
t84 = sin(qJ(2));
t138 = t81 * t84;
t80 = qJ(3) + qJ(4);
t74 = sin(t80);
t75 = cos(t80);
t161 = t124 * t75 - t74 * t138;
t85 = sin(qJ(1));
t137 = t81 * t85;
t114 = t85 * t124;
t88 = cos(qJ(2));
t89 = cos(qJ(1));
t54 = -t114 * t84 + t88 * t89;
t23 = t75 * t137 - t54 * t74;
t87 = cos(qJ(3));
t77 = t87 * pkin(3);
t73 = t77 + pkin(2);
t83 = sin(qJ(3));
t153 = m(4) * pkin(2) + m(5) * t73 + t87 * mrSges(4,1) + t75 * mrSges(5,1) - t83 * mrSges(4,2) - t74 * mrSges(5,2) + mrSges(3,1);
t76 = qJ(5) + t80;
t71 = sin(t76);
t72 = cos(t76);
t148 = -t110 * t71 + t164 * t72 + t153;
t158 = -m(5) * pkin(3) - mrSges(4,1);
t38 = t124 * t72 - t138 * t71;
t39 = t124 * t71 + t138 * t72;
t156 = -t162 * t38 + t163 * t39;
t17 = -t137 * t72 + t54 * t71;
t18 = t137 * t71 + t54 * t72;
t155 = t162 * t17 + t163 * t18;
t134 = t81 * t89;
t113 = t89 * t124;
t52 = t113 * t84 + t85 * t88;
t13 = -t72 * t134 - t52 * t71;
t14 = -t71 * t134 + t52 * t72;
t154 = -t162 * t13 + t163 * t14;
t151 = -t161 * mrSges(5,1) - (-t124 * t74 - t138 * t75) * mrSges(5,2) + t156;
t24 = t137 * t74 + t54 * t75;
t150 = -t23 * mrSges(5,1) + t24 * mrSges(5,2) + t155;
t100 = -t134 * t75 - t52 * t74;
t63 = t74 * t134;
t149 = -t100 * mrSges(5,1) - (-t52 * t75 + t63) * mrSges(5,2) + t154;
t144 = pkin(3) * t83;
t136 = t81 * t87;
t135 = t81 * t88;
t61 = pkin(4) * t74 + t144;
t62 = pkin(4) * t75 + t77;
t129 = t62 * t137 - t54 * t61;
t126 = t124 * t62 - t61 * t138;
t125 = t89 * pkin(1) + pkin(8) * t137;
t118 = t85 * pkin(1) - pkin(8) * t134;
t117 = t13 * pkin(5) + pkin(12) * t14;
t116 = -t17 * pkin(5) + pkin(12) * t18;
t115 = t38 * pkin(5) + pkin(12) * t39;
t111 = -m(5) * t144 - mrSges(3,3);
t109 = t23 * pkin(4);
t108 = -t134 * t62 - t52 * t61;
t53 = t114 * t88 + t89 * t84;
t59 = pkin(2) + t62;
t79 = -pkin(11) + t90;
t107 = t61 * t137 - t53 * t79 + t54 * t59 + t125;
t102 = t161 * pkin(4);
t25 = t136 * t85 - t54 * t83;
t95 = t100 * pkin(4);
t66 = t83 * t134;
t51 = -t113 * t88 + t84 * t85;
t26 = t137 * t83 + t54 * t87;
t2 = t18 * t86 + t53 * t82;
t1 = -t18 * t82 + t53 * t86;
t3 = [(-t89 * mrSges(2,1) - m(3) * t125 - t54 * mrSges(3,1) - m(4) * (pkin(2) * t54 + t125) - t26 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * (t54 * t73 + t125) - t24 * mrSges(5,1) - t23 * mrSges(5,2) - m(6) * t107 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t107) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (t111 * t81 + mrSges(2,2)) * t85 + t110 * t17 + t92 * t53) * g(2) + (t85 * mrSges(2,1) - t66 * mrSges(4,1) - t63 * mrSges(5,1) + t153 * t52 + (mrSges(2,2) + (-mrSges(4,2) * t87 - mrSges(5,2) * t75 + t111) * t81) * t89 + t110 * t13 + t164 * t14 + t157 * t51 + (m(3) + m(4) + m(5)) * t118 + t159 * (-t61 * t134 - t51 * t79 + t52 * t59 + t118)) * g(1) (-t159 * (-t51 * t59 - t52 * t79) - t157 * t52 + t148 * t51) * g(2) + (-t159 * (-t53 * t59 - t54 * t79) - t157 * t54 + t148 * t53) * g(1) + (-t159 * t59 * t135 + (-t148 * t88 + (t159 * t79 - t157) * t84) * t81) * g(3) (-(-t124 * t83 - t136 * t84) * mrSges(4,2) - m(6) * t126 - m(7) * (t115 + t126) + t158 * (t124 * t87 - t138 * t83) + t151) * g(3) + (-(-t52 * t87 + t66) * mrSges(4,2) - m(6) * t108 - m(7) * (t108 + t117) + t158 * (-t134 * t87 - t52 * t83) + t149) * g(2) + (mrSges(4,2) * t26 - m(6) * t129 - m(7) * (t116 + t129) + t158 * t25 + t150) * g(1) (-m(6) * t102 - m(7) * (t102 + t115) + t151) * g(3) + (-m(6) * t95 - m(7) * (t117 + t95) + t149) * g(2) + (-m(6) * t109 - m(7) * (t109 + t116) + t150) * g(1) (-m(7) * t115 + t156) * g(3) + (-m(7) * t117 + t154) * g(2) + (-m(7) * t116 + t155) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t82 + t51 * t86) * mrSges(7,1) + (-t14 * t86 - t51 * t82) * mrSges(7,2)) - g(3) * ((-t135 * t86 - t39 * t82) * mrSges(7,1) + (t135 * t82 - t39 * t86) * mrSges(7,2))];
taug  = t3(:);
