% Calculate Gravitation load on the joints for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:36
% EndTime: 2019-03-10 01:34:41
% DurationCPUTime: 1.81s
% Computational Cost: add. (954->164), mult. (1619->218), div. (0->0), fcn. (1882->12), ass. (0->78)
t173 = mrSges(6,2) + mrSges(7,2);
t183 = mrSges(6,1) + mrSges(7,1);
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t190 = -t173 * t85 + t89 * t183 + mrSges(5,1);
t181 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t84 = -qJ(6) - pkin(11);
t189 = m(6) * pkin(11) - m(7) * t84 - t181;
t169 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t159 = cos(qJ(1));
t83 = sin(pkin(6));
t116 = t83 * t159;
t121 = cos(pkin(6));
t104 = t121 * t159;
t87 = sin(qJ(2));
t88 = sin(qJ(1));
t91 = cos(qJ(2));
t61 = t87 * t104 + t88 * t91;
t82 = qJ(3) + qJ(4);
t79 = sin(t82);
t80 = cos(t82);
t31 = t80 * t116 + t61 * t79;
t32 = -t79 * t116 + t61 * t80;
t187 = t181 * t32 + t190 * t31;
t142 = t83 * t88;
t109 = t88 * t121;
t63 = -t87 * t109 + t159 * t91;
t35 = -t80 * t142 + t63 * t79;
t36 = t79 * t142 + t63 * t80;
t186 = t181 * t36 + t190 * t35;
t143 = t83 * t87;
t54 = -t121 * t80 + t79 * t143;
t55 = t121 * t79 + t80 * t143;
t185 = t181 * t55 + t190 * t54;
t162 = m(7) * pkin(5);
t96 = t85 * t162 - t169;
t77 = pkin(5) * t89 + pkin(4);
t86 = sin(qJ(3));
t90 = cos(qJ(3));
t163 = m(4) * pkin(2) + t90 * mrSges(4,1) - t86 * mrSges(4,2) + mrSges(3,1) + (m(6) * pkin(4) + m(7) * t77 + mrSges(5,1)) * t80 + t189 * t79;
t178 = t121 * t90 - t86 * t143;
t141 = t83 * t90;
t37 = t88 * t141 - t63 * t86;
t60 = -t91 * t104 + t87 * t88;
t176 = t32 * t85 - t60 * t89;
t175 = -t32 * t89 - t60 * t85;
t171 = m(5) + m(6) + m(7);
t170 = -t162 - t183;
t150 = t61 * t85;
t149 = t63 * t85;
t145 = t80 * t85;
t144 = t80 * t89;
t140 = t85 * t91;
t139 = t89 * t91;
t136 = -t31 * t77 - t32 * t84;
t135 = -t35 * t77 - t36 * t84;
t128 = -t54 * t77 - t55 * t84;
t122 = t159 * pkin(1) + pkin(8) * t142;
t120 = t86 * t142;
t115 = -pkin(1) * t88 + pkin(8) * t116;
t114 = -t31 * pkin(4) + t32 * pkin(11);
t113 = -t35 * pkin(4) + pkin(11) * t36;
t112 = -t54 * pkin(4) + pkin(11) * t55;
t71 = t86 * t116;
t110 = -t61 * t90 + t71;
t107 = t37 * pkin(3);
t62 = t91 * t109 + t159 * t87;
t78 = pkin(3) * t90 + pkin(2);
t92 = -pkin(10) - pkin(9);
t106 = pkin(3) * t120 - t62 * t92 + t63 * t78 + t122;
t5 = -t36 * t85 + t62 * t89;
t101 = t178 * pkin(3);
t100 = pkin(3) * t71 + t60 * t92 - t61 * t78 + t115;
t98 = t90 * t116 + t61 * t86;
t95 = t98 * pkin(3);
t38 = t63 * t90 + t120;
t6 = t36 * t89 + t62 * t85;
t1 = [(-t159 * mrSges(2,1) - m(3) * t122 - t63 * mrSges(3,1) - m(4) * (pkin(2) * t63 + t122) - t38 * mrSges(4,1) - t37 * mrSges(4,2) - m(5) * t106 - t36 * mrSges(5,1) - m(6) * (pkin(4) * t36 + t106) - m(7) * (t36 * t77 + t106) + (-mrSges(3,3) * t83 + mrSges(2,2)) * t88 - t183 * t6 - t173 * t5 - t96 * t62 - t189 * t35) * g(2) + (t88 * mrSges(2,1) + t159 * mrSges(2,2) - m(3) * t115 + t61 * mrSges(3,1) - mrSges(3,3) * t116 - m(4) * (-pkin(2) * t61 + t115) - t110 * mrSges(4,1) - t98 * mrSges(4,2) - m(5) * t100 + t32 * mrSges(5,1) - m(6) * (-pkin(4) * t32 + t100) - m(7) * (-t32 * t77 + t100) - t183 * t175 - t173 * t176 + t96 * t60 + t189 * t31) * g(1) (-t150 * t162 - t183 * (-t60 * t144 + t150) - t173 * (t60 * t145 + t61 * t89) - t171 * (-t60 * t78 - t61 * t92) + t169 * t61 + t163 * t60) * g(2) + (-t149 * t162 - t173 * (t62 * t145 + t63 * t89) - t171 * (-t62 * t78 - t63 * t92) - t183 * (-t62 * t144 + t149) + t169 * t63 + t163 * t62) * g(1) + ((-t171 * t78 - t163) * t91 + (-t139 * t183 + t173 * t140) * t80 + (t171 * t92 - t173 * t89 - t183 * t85 - t96) * t87) * g(3) * t83 (-t178 * mrSges(4,1) - (-t121 * t86 - t87 * t141) * mrSges(4,2) - m(5) * t101 - m(6) * (t101 + t112) - m(7) * (t101 + t128) + t185) * g(3) + (t98 * mrSges(4,1) - t110 * mrSges(4,2) + m(5) * t95 - m(6) * (t114 - t95) - m(7) * (-t95 + t136) + t187) * g(2) + (-mrSges(4,1) * t37 + mrSges(4,2) * t38 - m(5) * t107 - m(6) * (t107 + t113) - m(7) * (t107 + t135) + t186) * g(1) (-m(6) * t112 - m(7) * t128 + t185) * g(3) + (-m(6) * t114 - m(7) * t136 + t187) * g(2) + (-m(6) * t113 - m(7) * t135 + t186) * g(1) (-t173 * (t83 * t140 - t55 * t89) + t170 * (-t83 * t139 - t55 * t85)) * g(3) + (-t170 * t176 - t173 * t175) * g(2) + (t170 * t5 + t173 * t6) * g(1) (-g(1) * t35 - g(2) * t31 - g(3) * t54) * m(7)];
taug  = t1(:);
