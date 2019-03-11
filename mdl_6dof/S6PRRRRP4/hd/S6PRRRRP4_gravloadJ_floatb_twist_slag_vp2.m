% Calculate Gravitation load on the joints for
% S6PRRRRP4
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:09
% EndTime: 2019-03-09 00:13:12
% DurationCPUTime: 1.08s
% Computational Cost: add. (729->120), mult. (1479->166), div. (0->0), fcn. (1761->12), ass. (0->75)
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t157 = -m(5) * pkin(3) - t74 * mrSges(5,1) + t71 * mrSges(5,2) - mrSges(4,1);
t142 = -m(6) - m(7);
t65 = pkin(4) * t74 + pkin(3);
t156 = t142 * t65;
t155 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(5,3);
t152 = mrSges(6,2) - mrSges(7,3);
t131 = -m(7) * qJ(6) + t152;
t153 = mrSges(6,1) + mrSges(7,1);
t132 = -m(7) * pkin(5) - t153;
t68 = qJ(4) + qJ(5);
t66 = sin(t68);
t67 = cos(t68);
t154 = -t131 * t66 - t132 * t67;
t151 = mrSges(6,3) + mrSges(7,2);
t150 = -m(4) + t142;
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t149 = t155 * t72 + t157 * t75 - mrSges(3,1);
t148 = -t71 * mrSges(5,1) - t74 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t103 = cos(pkin(11));
t73 = sin(qJ(2));
t76 = cos(qJ(2));
t104 = cos(pkin(6));
t69 = sin(pkin(11));
t96 = t69 * t104;
t53 = t103 * t76 - t73 * t96;
t70 = sin(pkin(6));
t33 = t69 * t70 * t72 + t53 * t75;
t52 = t103 * t73 + t76 * t96;
t147 = -t33 * t71 + t52 * t74;
t87 = t104 * t103;
t51 = t69 * t76 + t73 * t87;
t95 = t70 * t103;
t31 = t51 * t75 - t72 * t95;
t50 = t69 * t73 - t76 * t87;
t146 = -t31 * t71 + t50 * t74;
t125 = pkin(4) * t71;
t145 = -m(5) * pkin(8) + t142 * t125 - t131 * t67 + t132 * t66 + t148;
t77 = -pkin(10) - pkin(9);
t144 = -t149 + (t142 * t77 + t151) * t72 + (t154 - t156) * t75;
t143 = -m(4) - m(5);
t115 = t70 * t76;
t116 = t70 * t75;
t55 = t104 * t72 + t116 * t73;
t26 = t115 * t67 + t55 * t66;
t101 = t66 * t115;
t27 = t55 * t67 - t101;
t138 = t152 * t27 + t153 * t26;
t13 = t33 * t66 - t52 * t67;
t14 = t33 * t67 + t52 * t66;
t137 = t153 * t13 + t152 * t14;
t11 = t31 * t66 - t50 * t67;
t12 = t31 * t67 + t50 * t66;
t136 = t153 * t11 + t152 * t12;
t134 = -t154 + t157;
t133 = -t151 + t155;
t117 = t70 * t73;
t111 = t75 * t76;
t105 = pkin(2) * t115 + pkin(8) * t117;
t102 = t72 * t115;
t97 = -t11 * pkin(5) + qJ(6) * t12;
t94 = -t13 * pkin(5) + qJ(6) * t14;
t93 = -t26 * pkin(5) + qJ(6) * t27;
t92 = t146 * pkin(4);
t91 = t147 * pkin(4);
t86 = -t115 * t74 - t55 * t71;
t83 = t86 * pkin(4);
t54 = t104 * t75 - t117 * t72;
t49 = t52 * pkin(2);
t48 = t50 * pkin(2);
t32 = t116 * t69 - t53 * t72;
t30 = -t51 * t72 - t75 * t95;
t1 = [(-m(2) - m(3) + t142 + t143) * g(3) (t142 * (-t102 * t77 + t117 * t125 + t105) + t131 * (t101 * t75 - t117 * t67) + t143 * t105 - t151 * t102 + (t111 * t156 + t132 * (t111 * t67 + t66 * t73) + t149 * t76 + t148 * t73) * t70) * g(3) + (m(5) * t48 + t150 * (t51 * pkin(8) - t48) + t145 * t51 + t144 * t50) * g(2) + (m(5) * t49 + t150 * (t53 * pkin(8) - t49) + t145 * t53 + t144 * t52) * g(1) (t142 * (t54 * t65 - t55 * t77) + t133 * t55 + t134 * t54) * g(3) + (t142 * (t30 * t65 - t31 * t77) + t133 * t31 + t134 * t30) * g(2) + (t142 * (t32 * t65 - t33 * t77) + t133 * t33 + t134 * t32) * g(1) (-t86 * mrSges(5,1) - (t115 * t71 - t55 * t74) * mrSges(5,2) - m(6) * t83 - m(7) * (t83 + t93) + t138) * g(3) + (-t146 * mrSges(5,1) - (-t31 * t74 - t50 * t71) * mrSges(5,2) - m(6) * t92 - m(7) * (t92 + t97) + t136) * g(2) + (-t147 * mrSges(5,1) - (-t33 * t74 - t52 * t71) * mrSges(5,2) - m(6) * t91 - m(7) * (t91 + t94) + t137) * g(1) (-m(7) * t93 + t138) * g(3) + (-m(7) * t97 + t136) * g(2) + (-m(7) * t94 + t137) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t26) * m(7)];
taug  = t1(:);
