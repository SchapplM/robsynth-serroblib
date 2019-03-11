% Calculate Gravitation load on the joints for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:45
% EndTime: 2019-03-10 05:16:51
% DurationCPUTime: 2.13s
% Computational Cost: add. (1693->177), mult. (4412->256), div. (0->0), fcn. (5605->16), ass. (0->90)
t149 = cos(pkin(6));
t155 = cos(qJ(2));
t134 = t149 * t155;
t152 = sin(qJ(2));
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t106 = -t134 * t156 + t153 * t152;
t146 = sin(pkin(7));
t147 = sin(pkin(6));
t121 = t147 * t146;
t148 = cos(pkin(7));
t181 = t106 * t148 + t156 * t121;
t154 = cos(qJ(3));
t133 = t149 * t152;
t60 = t133 * t156 + t153 * t155;
t82 = sin(qJ(3));
t30 = -t154 * t60 + t181 * t82;
t122 = t148 * t147;
t46 = t106 * t146 - t156 * t122;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t180 = t30 * t84 - t46 * t81;
t179 = t30 * t81 + t46 * t84;
t176 = -m(5) - m(6);
t80 = sin(qJ(5));
t178 = mrSges(4,2) - mrSges(5,3) - m(7) * (pkin(5) * t80 + pkin(11));
t83 = cos(qJ(5));
t75 = pkin(5) * t83 + pkin(4);
t79 = qJ(5) + qJ(6);
t76 = sin(t79);
t77 = cos(t79);
t168 = m(6) * pkin(4) + m(7) * t75 + t83 * mrSges(6,1) + t77 * mrSges(7,1) - t80 * mrSges(6,2) - t76 * mrSges(7,2) + mrSges(5,1);
t163 = mrSges(5,2) + m(7) * (-pkin(13) - pkin(12)) - mrSges(7,3) - m(6) * pkin(12) - mrSges(6,3);
t98 = t134 * t153 + t152 * t156;
t86 = t153 * t122 + t98 * t146;
t177 = mrSges(4,1) + t168 * t84 - t163 * t81 + pkin(3) * (m(7) - t176);
t162 = -t80 * mrSges(6,1) - t76 * mrSges(7,1) - t83 * mrSges(6,2) - t77 * mrSges(7,2) + t178;
t27 = t181 * t154 + t60 * t82;
t173 = -m(7) * pkin(5) - mrSges(6,1);
t172 = -t146 * mrSges(4,3) + mrSges(3,2);
t171 = t153 * t121 - t98 * t148;
t169 = t155 * t122 + t149 * t146;
t164 = t176 * pkin(11) + t162;
t159 = (t180 * t76 + t27 * t77) * mrSges(7,1) + (t180 * t77 - t27 * t76) * mrSges(7,2);
t61 = -t133 * t153 + t155 * t156;
t32 = t61 * t154 + t171 * t82;
t16 = t32 * t84 + t81 * t86;
t31 = -t171 * t154 + t61 * t82;
t5 = -t16 * t76 + t31 * t77;
t6 = t16 * t77 + t31 * t76;
t158 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t128 = t147 * t152;
t45 = t154 * t128 + t169 * t82;
t59 = -t121 * t155 + t148 * t149;
t26 = t45 * t84 + t59 * t81;
t44 = t128 * t82 - t169 * t154;
t157 = (-t26 * t76 + t44 * t77) * mrSges(7,1) + (-t26 * t77 - t44 * t76) * mrSges(7,2);
t109 = t152 * t121;
t130 = t155 * t147;
t151 = pkin(2) * t130 + pkin(10) * t109;
t129 = t147 * t153;
t150 = t156 * pkin(1) + pkin(9) * t129;
t113 = t152 * t122;
t54 = -t113 * t82 + t130 * t154;
t145 = t54 * pkin(3) + t151;
t141 = pkin(10) * t146;
t140 = t81 * t146;
t139 = t82 * t148;
t138 = t84 * t146;
t131 = t156 * t147;
t136 = -pkin(1) * t153 + pkin(9) * t131;
t132 = t148 * t154;
t7 = -t16 * t80 + t31 * t83;
t120 = -t106 * pkin(2) + t141 * t60;
t119 = -t98 * pkin(2) + t141 * t61;
t36 = -t106 * t154 - t139 * t60;
t118 = t36 * pkin(3) + t120;
t38 = -t139 * t61 - t154 * t98;
t117 = t38 * pkin(3) + t119;
t93 = -t60 * pkin(2) - pkin(10) * t46 + t136;
t91 = t30 * pkin(3) + t93;
t90 = t61 * pkin(2) + pkin(10) * t86 + t150;
t89 = t32 * pkin(3) + t90;
t87 = t31 * pkin(11) + t89;
t53 = t113 * t154 + t130 * t82;
t37 = t132 * t61 - t82 * t98;
t35 = -t106 * t82 + t132 * t60;
t15 = t32 * t81 - t84 * t86;
t8 = t16 * t83 + t31 * t80;
t1 = [(-t156 * mrSges(2,1) + t153 * mrSges(2,2) - m(3) * t150 - t61 * mrSges(3,1) + t98 * mrSges(3,2) - mrSges(3,3) * t129 - m(4) * t90 - t32 * mrSges(4,1) - t86 * mrSges(4,3) - m(5) * t87 - t16 * mrSges(5,1) - m(6) * (t16 * pkin(4) + t87) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t75 + t89) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t178 * t31 + t163 * t15) * g(2) + (-m(3) * t136 - m(4) * t93 - m(7) * t91 + t153 * mrSges(2,1) + t60 * mrSges(3,1) - t30 * mrSges(4,1) + t156 * mrSges(2,2) - t106 * mrSges(3,2) - mrSges(3,3) * t131 + t46 * mrSges(4,3) + t176 * (-pkin(11) * t27 + t91) - t168 * t180 - t162 * t27 + t163 * t179) * g(1) (-m(4) * t151 - m(7) * t145 - mrSges(3,1) * t130 - t54 * mrSges(4,1) + mrSges(3,2) * t128 - mrSges(4,3) * t109 + t176 * (t53 * pkin(11) + t145) - t168 * (t109 * t81 + t54 * t84) + t162 * t53 + t163 * (-t109 * t84 + t54 * t81)) * g(3) + (-m(4) * t120 - m(7) * t118 + mrSges(3,1) * t106 - t36 * mrSges(4,1) + t172 * t60 + t176 * (t35 * pkin(11) + t118) - t168 * (t140 * t60 + t36 * t84) + t162 * t35 + t163 * (-t138 * t60 + t36 * t81)) * g(2) + (-m(4) * t119 - m(7) * t117 + mrSges(3,1) * t98 - t38 * mrSges(4,1) + t172 * t61 + t176 * (t37 * pkin(11) + t117) - t168 * (t140 * t61 + t38 * t84) + t162 * t37 + t163 * (-t138 * t61 + t38 * t81)) * g(1) (t164 * t45 + t177 * t44) * g(3) + (-t164 * t30 + t177 * t27) * g(2) + (t164 * t32 + t177 * t31) * g(1) (t163 * t26 - t168 * (-t45 * t81 + t59 * t84)) * g(3) + (-t163 * t180 - t168 * t179) * g(2) + (t168 * t15 + t163 * t16) * g(1) (-(-t26 * t83 - t44 * t80) * mrSges(6,2) - t157 + t173 * (-t26 * t80 + t44 * t83)) * g(3) + (-(t180 * t83 - t27 * t80) * mrSges(6,2) - t159 + t173 * (t180 * t80 + t27 * t83)) * g(2) + (mrSges(6,2) * t8 + t173 * t7 - t158) * g(1), -g(1) * t158 - g(2) * t159 - g(3) * t157];
taug  = t1(:);
