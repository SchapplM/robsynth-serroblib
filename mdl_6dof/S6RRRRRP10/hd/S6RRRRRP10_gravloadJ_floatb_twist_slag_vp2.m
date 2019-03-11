% Calculate Gravitation load on the joints for
% S6RRRRRP10
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:20
% EndTime: 2019-03-10 02:16:25
% DurationCPUTime: 1.73s
% Computational Cost: add. (979->157), mult. (1991->207), div. (0->0), fcn. (2391->12), ass. (0->85)
t189 = mrSges(6,2) - mrSges(7,3);
t118 = -m(7) * qJ(6) + t189;
t190 = mrSges(6,1) + mrSges(7,1);
t195 = m(7) * pkin(5) + t190;
t88 = qJ(4) + qJ(5);
t85 = sin(t88);
t86 = cos(t88);
t197 = -t118 * t85 + t195 * t86;
t179 = m(6) + m(7);
t93 = cos(qJ(4));
t84 = pkin(4) * t93 + pkin(3);
t192 = t179 * t84;
t90 = sin(qJ(4));
t176 = -m(5) * pkin(3) - t93 * mrSges(5,1) + t90 * mrSges(5,2) - mrSges(4,1);
t191 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3);
t188 = mrSges(6,3) + mrSges(7,2);
t187 = -m(4) - t179;
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t186 = t176 * t94 + t191 * t91 - mrSges(3,1);
t185 = -t176 + t197;
t113 = t90 * mrSges(5,1) + t93 * mrSges(5,2);
t143 = mrSges(3,2) - mrSges(4,3);
t184 = -t113 + t143;
t161 = sin(qJ(1));
t89 = sin(pkin(6));
t128 = t89 * t161;
t135 = cos(pkin(6));
t115 = t135 * t161;
t162 = cos(qJ(1));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t70 = -t115 * t92 + t162 * t95;
t43 = t128 * t91 + t70 * t94;
t69 = t115 * t95 + t162 * t92;
t17 = -t43 * t90 + t69 * t93;
t129 = t89 * t162;
t116 = t135 * t162;
t68 = t116 * t92 + t161 * t95;
t39 = -t91 * t129 + t68 * t94;
t67 = -t116 * t95 + t161 * t92;
t183 = -t39 * t90 + t67 * t93;
t11 = t39 * t85 - t67 * t86;
t12 = t39 * t86 + t67 * t85;
t163 = pkin(4) * t90;
t182 = -m(5) * pkin(9) - t118 * t86 - t179 * t163 - t195 * t85 + t184;
t96 = -pkin(11) - pkin(10);
t181 = -t186 + (-t179 * t96 + t188) * t91 + (t192 + t197) * t94;
t180 = m(4) + m(5);
t148 = t89 * t95;
t149 = t89 * t92;
t66 = t135 * t91 + t149 * t94;
t32 = t148 * t86 + t66 * t85;
t133 = t85 * t148;
t33 = t66 * t86 - t133;
t175 = t189 * t33 + t190 * t32;
t15 = t43 * t85 - t69 * t86;
t16 = t43 * t86 + t69 * t85;
t174 = t190 * t15 + t189 * t16;
t173 = t190 * t11 + t189 * t12;
t170 = -pkin(9) * (t179 + t180) + t143;
t106 = -t188 + t191;
t156 = t67 * t90;
t154 = t69 * t90;
t144 = t94 * t95;
t137 = pkin(2) * t148 + pkin(9) * t149;
t136 = t162 * pkin(1) + pkin(8) * t128;
t134 = t91 * t148;
t132 = t70 * pkin(2) + t136;
t125 = -t11 * pkin(5) + qJ(6) * t12;
t124 = -t94 * t129 - t68 * t91;
t123 = -t15 * pkin(5) + qJ(6) * t16;
t122 = -t32 * pkin(5) + qJ(6) * t33;
t120 = t183 * pkin(4);
t119 = t17 * pkin(4);
t117 = -pkin(1) * t161 + pkin(8) * t129;
t109 = t68 * pkin(2) - t117;
t107 = -t148 * t93 - t66 * t90;
t103 = t107 * pkin(4);
t65 = t135 * t94 - t149 * t91;
t63 = t69 * pkin(2);
t61 = t67 * pkin(2);
t42 = -t128 * t94 + t70 * t91;
t18 = t43 * t93 + t154;
t1 = [(-t162 * mrSges(2,1) + t161 * mrSges(2,2) - m(3) * t136 - t70 * mrSges(3,1) - mrSges(3,3) * t128 - m(4) * t132 - t43 * mrSges(4,1) - m(5) * (pkin(3) * t43 + t132) - t18 * mrSges(5,1) - t17 * mrSges(5,2) + t170 * t69 - t195 * t16 + t118 * t15 + t106 * t42 - t179 * (pkin(4) * t154 - t42 * t96 + t43 * t84 + t132)) * g(2) + (t161 * mrSges(2,1) + t162 * mrSges(2,2) - m(3) * t117 + t68 * mrSges(3,1) - mrSges(3,3) * t129 - t176 * t39 + (t113 - t170) * t67 + t195 * t12 - t118 * t11 + t106 * t124 + t180 * t109 + t179 * (pkin(4) * t156 + t124 * t96 + t39 * t84 + t109)) * g(1) (-t179 * (-t134 * t96 + t149 * t163 + t137) + t118 * (t133 * t94 - t149 * t86) - t180 * t137 - t188 * t134 + (-t144 * t192 - t195 * (t144 * t86 + t85 * t92) + t186 * t95 + t184 * t92) * t89) * g(3) + (m(5) * t61 + t187 * (t68 * pkin(9) - t61) + t182 * t68 + t181 * t67) * g(2) + (m(5) * t63 + t187 * (t70 * pkin(9) - t63) + t182 * t70 + t181 * t69) * g(1) (-t179 * (t65 * t84 - t66 * t96) + t106 * t66 - t185 * t65) * g(3) + (-t179 * (t124 * t84 - t39 * t96) + t106 * t39 - t185 * t124) * g(2) + (-t179 * (-t42 * t84 - t43 * t96) + t106 * t43 + t185 * t42) * g(1) (-t107 * mrSges(5,1) - (t148 * t90 - t66 * t93) * mrSges(5,2) - m(6) * t103 - m(7) * (t103 + t122) + t175) * g(3) + (-t183 * mrSges(5,1) - (-t39 * t93 - t156) * mrSges(5,2) - m(6) * t120 - m(7) * (t120 + t125) + t173) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t119 - m(7) * (t119 + t123) + t174) * g(1) (-m(7) * t122 + t175) * g(3) + (-m(7) * t125 + t173) * g(2) + (-m(7) * t123 + t174) * g(1) (-g(1) * t15 - g(2) * t11 - g(3) * t32) * m(7)];
taug  = t1(:);
