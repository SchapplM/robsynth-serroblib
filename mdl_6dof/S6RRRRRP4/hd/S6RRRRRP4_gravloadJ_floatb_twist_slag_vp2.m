% Calculate Gravitation load on the joints for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:07
% EndTime: 2019-03-10 01:12:09
% DurationCPUTime: 1.15s
% Computational Cost: add. (677->148), mult. (715->167), div. (0->0), fcn. (664->10), ass. (0->84)
t171 = mrSges(6,3) + mrSges(7,2);
t170 = mrSges(6,1) + mrSges(7,1);
t169 = mrSges(6,2) - mrSges(7,3);
t63 = sin(qJ(4));
t120 = t63 * mrSges(5,2);
t61 = qJ(4) + qJ(5);
t56 = sin(t61);
t62 = qJ(2) + qJ(3);
t57 = sin(t62);
t130 = t56 * t57;
t167 = -mrSges(6,2) * t130 - t120 * t57;
t166 = -mrSges(5,3) - t171;
t165 = m(6) + m(7);
t59 = cos(t62);
t155 = pkin(3) * t59 + pkin(9) * t57;
t164 = m(5) * t155;
t163 = t171 * t57;
t162 = -m(7) * qJ(6) - mrSges(7,3);
t161 = -t59 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t57;
t65 = sin(qJ(1));
t66 = cos(qJ(4));
t116 = t65 * t66;
t68 = cos(qJ(1));
t118 = t63 * t68;
t17 = -t118 * t59 + t116;
t160 = g(1) * t68 + g(2) * t65;
t106 = qJ(6) * t56;
t115 = t66 * mrSges(5,1);
t58 = cos(t61);
t127 = t57 * t58;
t54 = pkin(4) * t66 + pkin(3);
t159 = mrSges(6,1) * t127 + (t115 - m(7) * (-pkin(5) * t58 - t106 - t54) + t58 * mrSges(7,1) + t56 * mrSges(7,3)) * t57;
t144 = m(7) * pkin(5) + t170;
t158 = mrSges(6,2) * t58 + t144 * t56;
t113 = t68 * t56;
t117 = t65 * t58;
t13 = t113 * t59 - t117;
t122 = t59 * t68;
t14 = t122 * t58 + t56 * t65;
t152 = t13 * t170 + t14 * t169;
t123 = t59 * t65;
t11 = t123 * t56 + t58 * t68;
t12 = t117 * t59 - t113;
t151 = t11 * t170 + t12 * t169;
t150 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t69 = -pkin(10) - pkin(9);
t121 = t59 * t69;
t136 = pkin(3) * t57;
t64 = sin(qJ(2));
t137 = pkin(2) * t64;
t85 = -t54 * t57 - t121;
t149 = -m(7) * (-t121 - t137) - m(6) * (t85 - t137) - m(5) * (-t136 - t137) + t159;
t148 = -m(6) * t85 + m(7) * t121 + t159;
t147 = t122 * t166 + t167 * t68;
t67 = cos(qJ(2));
t89 = t67 * mrSges(3,1) - t64 * mrSges(3,2);
t146 = m(3) * pkin(1) + mrSges(2,1) - t161 + t89;
t124 = t58 * t59;
t145 = -t170 * t124 + t161 - t163 + (t169 * t56 - t115 + t120) * t59;
t143 = -mrSges(6,2) - t162;
t142 = t167 * t65 + (-m(5) * pkin(9) + t166) * t123;
t60 = t67 * pkin(2);
t126 = t57 * t68;
t125 = t57 * t69;
t26 = t59 * t54;
t119 = t63 * t65;
t114 = t66 * t68;
t70 = -pkin(8) - pkin(7);
t112 = t68 * t70;
t98 = t162 * t127;
t96 = -pkin(5) * t11 + qJ(6) * t12;
t95 = t26 - t125;
t55 = t60 + pkin(1);
t94 = t55 * t68 - t65 * t70;
t92 = -pkin(5) * t13 + qJ(6) * t14;
t86 = mrSges(4,1) * t57 + mrSges(4,2) * t59;
t84 = t17 * pkin(4);
t81 = pkin(5) * t124 + t106 * t59 + t95;
t15 = t119 * t59 + t114;
t79 = t15 * pkin(4);
t42 = pkin(9) * t122;
t18 = t114 * t59 + t119;
t16 = -t116 * t59 + t118;
t1 = [(-t18 * mrSges(5,1) - t17 * mrSges(5,2) + (-m(4) - m(5)) * t94 - t165 * (pkin(4) * t119 + t54 * t122 - t68 * t125 + t94) - t144 * t14 - t143 * t13 - t171 * t126 + t150 * t65 + (-t146 - t164) * t68) * g(2) + (m(5) * t112 - t16 * mrSges(5,1) - t15 * mrSges(5,2) - t165 * (pkin(4) * t118 + t65 * t125 - t112) + t144 * t12 + t143 * t11 + (m(4) * t70 + t150) * t68 + (m(4) * t55 - m(5) * (-t55 - t155) - t165 * (-t55 - t26) + t146 + t163) * t65) * g(1) (t149 * t65 + t142) * g(2) + (-m(5) * t42 + t149 * t68 + t147) * g(1) + (-t89 - m(4) * t60 - m(5) * (t60 + t155) - m(6) * (t60 + t95) - m(7) * (t60 + t81) + t145) * g(3) + t160 * (m(4) * t137 + mrSges(3,1) * t64 + mrSges(3,2) * t67 + t86) t160 * t86 + ((m(5) * t136 + t148) * t65 + t142) * g(2) + (-m(5) * (-pkin(3) * t126 + t42) + t148 * t68 + t147) * g(1) + (-m(6) * t95 - m(7) * t81 + t145 - t164) * g(3) (t15 * mrSges(5,1) - t16 * mrSges(5,2) + m(6) * t79 - m(7) * (-t79 + t96) + t151) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t84 - m(7) * (t84 + t92) + t152) * g(1) + (t98 + (mrSges(5,2) * t66 + t158 + (pkin(4) * t165 + mrSges(5,1)) * t63) * t57) * g(3) (t158 * t57 + t98) * g(3) + (-m(7) * t96 + t151) * g(2) + (-m(7) * t92 + t152) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t130) * m(7)];
taug  = t1(:);
