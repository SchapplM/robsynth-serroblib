% Calculate Gravitation load on the joints for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:11
% EndTime: 2019-03-09 00:01:14
% DurationCPUTime: 1.02s
% Computational Cost: add. (795->124), mult. (1352->170), div. (0->0), fcn. (1583->12), ass. (0->77)
t170 = mrSges(6,2) - mrSges(7,3);
t81 = sin(qJ(5));
t172 = t81 * t170 - mrSges(5,1);
t171 = -mrSges(6,1) - mrSges(7,1);
t169 = mrSges(6,3) + mrSges(7,2);
t103 = m(7) * pkin(5) - t171;
t78 = qJ(3) + qJ(4);
t77 = cos(t78);
t168 = t103 * t77;
t167 = -m(4) * pkin(8) - t103 * t81 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t166 = mrSges(5,2) - t169;
t157 = -m(6) - m(7);
t164 = -m(5) + t157;
t121 = cos(pkin(6));
t80 = sin(pkin(6));
t83 = sin(qJ(2));
t131 = t80 * t83;
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t162 = t121 * t85 - t82 * t131;
t130 = t80 * t85;
t79 = sin(pkin(11));
t108 = t79 * t121;
t120 = cos(pkin(11));
t86 = cos(qJ(2));
t64 = -t83 * t108 + t120 * t86;
t161 = t79 * t130 - t64 * t82;
t76 = sin(t78);
t160 = -m(4) * pkin(2) - t85 * mrSges(4,1) - t77 * mrSges(5,1) + t82 * mrSges(4,2) + t76 * mrSges(5,2) - mrSges(3,1);
t101 = -m(7) * qJ(6) + t170;
t84 = cos(qJ(5));
t159 = -t101 * t84 + t167;
t146 = pkin(4) * t77;
t158 = t157 * (-pkin(10) * t76 - t146) + t84 * t168 - t101 * t77 * t81 - t160 + t169 * t76;
t122 = qJ(6) * t81;
t107 = t80 * t120;
t98 = t121 * t120;
t62 = t79 * t86 + t83 * t98;
t29 = -t77 * t107 - t62 * t76;
t142 = t29 * t84;
t156 = pkin(5) * t142 + t29 * t122;
t132 = t79 * t80;
t31 = t77 * t132 - t64 * t76;
t140 = t31 * t84;
t154 = pkin(5) * t140 + t31 * t122;
t53 = t121 * t77 - t76 * t131;
t138 = t53 * t84;
t153 = pkin(5) * t138 + t53 * t122;
t54 = t121 * t76 + t77 * t131;
t150 = t171 * t138 + t166 * t54 + t172 * t53;
t32 = t76 * t132 + t64 * t77;
t149 = t171 * t140 + t166 * t32 + t172 * t31;
t30 = -t76 * t107 + t62 * t77;
t148 = t171 * t142 + t166 * t30 + t172 * t29;
t129 = t80 * t86;
t128 = t84 * t86;
t119 = t76 * t129;
t116 = t81 * t129;
t111 = t29 * pkin(4) + t30 * pkin(10);
t110 = t31 * pkin(4) + pkin(10) * t32;
t109 = t53 * pkin(4) + pkin(10) * t54;
t102 = t161 * pkin(3);
t99 = t162 * pkin(3);
t92 = t102 + t110;
t91 = -t85 * t107 - t62 * t82;
t90 = t109 + t99;
t89 = t91 * pkin(3);
t88 = t111 + t89;
t87 = -pkin(9) - pkin(8);
t75 = pkin(3) * t85 + pkin(2);
t65 = t75 * t129;
t63 = t86 * t108 + t120 * t83;
t61 = t79 * t83 - t86 * t98;
t33 = t80 * t128 + t54 * t81;
t3 = t32 * t81 - t63 * t84;
t1 = t30 * t81 - t61 * t84;
t2 = [(-m(2) - m(3) - m(4) + t164) * g(3) (t164 * (-t61 * t75 - t62 * t87) + t159 * t62 + t158 * t61) * g(2) + (t164 * (-t63 * t75 - t64 * t87) + t159 * t64 + t158 * t63) * g(1) + (-m(5) * t65 + t157 * (pkin(10) * t119 + t129 * t146 - t87 * t131 + t65) + t101 * (t77 * t116 - t84 * t131) - t169 * t119 + (-t128 * t168 + t160 * t86 + (m(5) * t87 + t167) * t83) * t80) * g(3) (-t162 * mrSges(4,1) - (-t121 * t82 - t83 * t130) * mrSges(4,2) - m(5) * t99 - m(6) * t90 - m(7) * (t90 + t153) + t150) * g(3) + (-t91 * mrSges(4,1) - (t82 * t107 - t62 * t85) * mrSges(4,2) - m(5) * t89 - m(6) * t88 - m(7) * (t88 + t156) + t148) * g(2) + (-t161 * mrSges(4,1) - (-t82 * t132 - t64 * t85) * mrSges(4,2) - m(5) * t102 - m(6) * t92 - m(7) * (t92 + t154) + t149) * g(1) (-m(6) * t109 - m(7) * (t109 + t153) + t150) * g(3) + (-m(6) * t111 - m(7) * (t111 + t156) + t148) * g(2) + (-m(6) * t110 - m(7) * (t110 + t154) + t149) * g(1) (t101 * (t54 * t84 - t116) + t103 * t33) * g(3) + (t101 * (t30 * t84 + t61 * t81) + t103 * t1) * g(2) + (t101 * (t32 * t84 + t63 * t81) + t103 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t33) * m(7)];
taug  = t2(:);
