% Calculate Gravitation load on the joints for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:13
% EndTime: 2019-03-10 03:39:17
% DurationCPUTime: 1.08s
% Computational Cost: add. (683->143), mult. (664->165), div. (0->0), fcn. (603->12), ass. (0->84)
t63 = qJ(4) + qJ(5);
t59 = qJ(6) + t63;
t51 = sin(t59);
t55 = sin(t63);
t65 = sin(qJ(4));
t168 = -t65 * mrSges(5,2) - t55 * mrSges(6,2) - t51 * mrSges(7,2);
t167 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t52 = cos(t59);
t57 = cos(t63);
t68 = cos(qJ(4));
t166 = t68 * mrSges(5,1) + t57 * mrSges(6,1) + t52 * mrSges(7,1);
t64 = qJ(2) + qJ(3);
t56 = sin(t64);
t165 = t168 * t56;
t164 = -m(5) * pkin(9) + t167;
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t162 = g(1) * t70 + g(2) * t67;
t161 = t166 * t56;
t58 = cos(t64);
t160 = -t58 * mrSges(4,1) + (mrSges(4,2) + t167) * t56;
t60 = t68 * pkin(4);
t53 = t60 + pkin(3);
t71 = -pkin(10) - pkin(9);
t153 = t58 * t53 - t56 * t71;
t31 = pkin(5) * t57 + t60;
t27 = pkin(3) + t31;
t62 = -pkin(11) + t71;
t154 = t58 * t27 - t56 * t62;
t158 = t58 * pkin(3) + t56 * pkin(9);
t159 = -m(5) * t158 - m(6) * t153 - m(7) * t154;
t140 = m(6) * pkin(4);
t87 = -mrSges(7,1) * t51 - mrSges(7,2) * t52;
t157 = mrSges(6,1) * t55 + mrSges(6,2) * t57 - t87;
t155 = mrSges(5,1) + t140;
t120 = t58 * t70;
t7 = -t51 * t120 + t52 * t67;
t8 = t52 * t120 + t51 * t67;
t137 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t15 = -t55 * t120 + t57 * t67;
t16 = t57 * t120 + t55 * t67;
t151 = -t15 * mrSges(6,1) + t16 * mrSges(6,2) - t137;
t121 = t58 * t67;
t13 = t55 * t121 + t57 * t70;
t5 = t51 * t121 + t52 * t70;
t6 = -t52 * t121 + t51 * t70;
t138 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t14 = -t57 * t121 + t55 * t70;
t150 = t13 * mrSges(6,1) - t14 * mrSges(6,2) - t138;
t149 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t135 = pkin(3) * t56;
t66 = sin(qJ(2));
t136 = pkin(2) * t66;
t83 = -t53 * t56 - t58 * t71;
t85 = -t27 * t56 - t58 * t62;
t148 = -m(7) * (t85 - t136) - m(6) * (t83 - t136) - m(5) * (-t135 - t136) + t161;
t147 = m(4) + m(5) + m(6) + m(7);
t146 = t160 + (-t166 - t168) * t58;
t134 = pkin(4) * t65;
t133 = pkin(5) * t55;
t30 = t133 + t134;
t145 = -m(6) * t134 - m(7) * t30;
t144 = m(5) * t135 - m(6) * t83 - m(7) * t85 + t161;
t69 = cos(qJ(2));
t92 = t69 * mrSges(3,1) - t66 * mrSges(3,2);
t143 = m(3) * pkin(1) + mrSges(2,1) - t160 + t92;
t142 = t164 * t120 + t165 * t70;
t141 = t164 * t121 + t165 * t67;
t139 = m(7) * pkin(5);
t130 = g(3) * t56;
t61 = t69 * pkin(2);
t118 = t65 * t67;
t117 = t65 * t70;
t116 = t67 * t30;
t115 = t67 * t68;
t113 = t68 * t70;
t89 = mrSges(4,1) * t56 + mrSges(4,2) * t58;
t19 = -t58 * t117 + t115;
t17 = t58 * t118 + t113;
t72 = -pkin(8) - pkin(7);
t54 = t61 + pkin(1);
t20 = t58 * t113 + t118;
t18 = -t58 * t115 + t117;
t1 = [(-t118 * t140 - m(7) * t116 - t20 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t19 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) - t147 * (t70 * t54 - t67 * t72) + t149 * t67 + (-t143 + t159) * t70) * g(2) + (-t18 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t17 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (t147 * t72 + t145 + t149) * t70 + (m(4) * t54 - m(5) * (-t54 - t158) - m(6) * (-t54 - t153) - m(7) * (-t54 - t154) + t143) * t67) * g(1) (t148 * t67 + t141) * g(2) + (t148 * t70 + t142) * g(1) + (-t92 - m(4) * t61 - m(5) * (t61 + t158) - m(6) * (t61 + t153) - m(7) * (t61 + t154) + t146) * g(3) + t162 * (m(4) * t136 + mrSges(3,1) * t66 + mrSges(3,2) * t69 + t89) t162 * t89 + (t144 * t67 + t141) * g(2) + (t144 * t70 + t142) * g(1) + (t146 + t159) * g(3) (mrSges(5,1) * t65 + mrSges(5,2) * t68 - t145 + t157) * t130 + (-t18 * mrSges(5,2) - m(7) * (-t58 * t116 - t31 * t70) + t155 * t17 + t150) * g(2) + (t20 * mrSges(5,2) - m(7) * (-t30 * t120 + t31 * t67) - t155 * t19 + t151) * g(1) (m(7) * t133 + t157) * t130 + (t13 * t139 + t150) * g(2) + (-t15 * t139 + t151) * g(1), -g(1) * t137 - g(2) * t138 - t87 * t130];
taug  = t1(:);
