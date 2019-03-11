% Calculate Gravitation load on the joints for
% S6PRRRRP1
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:46
% EndTime: 2019-03-08 23:55:49
% DurationCPUTime: 1.10s
% Computational Cost: add. (712->121), mult. (1197->167), div. (0->0), fcn. (1378->12), ass. (0->63)
t155 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t140 = -mrSges(6,2) - mrSges(7,2);
t156 = mrSges(6,1) + mrSges(7,1);
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t157 = t140 * t69 + t72 * t156 + mrSges(5,1);
t138 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t132 = m(7) * pkin(5);
t137 = -t132 - t156;
t62 = pkin(5) * t72 + pkin(4);
t66 = qJ(3) + qJ(4);
t64 = sin(t66);
t65 = cos(t66);
t68 = -qJ(6) - pkin(10);
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t133 = m(4) * pkin(2) + t73 * mrSges(4,1) - t70 * mrSges(4,2) + mrSges(3,1) + (m(6) * pkin(4) + m(7) * t62 + mrSges(5,1)) * t65 + (m(6) * pkin(10) - m(7) * t68 - t155) * t64;
t100 = cos(pkin(6));
t67 = sin(pkin(6));
t71 = sin(qJ(2));
t116 = t67 * t71;
t152 = t100 * t73 - t70 * t116;
t74 = cos(qJ(2));
t98 = sin(pkin(11));
t81 = t100 * t98;
t99 = cos(pkin(11));
t53 = -t71 * t81 + t99 * t74;
t90 = t67 * t98;
t151 = -t53 * t70 + t73 * t90;
t82 = t100 * t99;
t51 = t71 * t82 + t98 * t74;
t91 = t67 * t99;
t27 = t51 * t64 + t65 * t91;
t28 = t51 * t65 - t64 * t91;
t150 = t155 * t28 + t157 * t27;
t29 = t53 * t64 - t65 * t90;
t30 = t53 * t65 + t64 * t90;
t149 = t155 * t30 + t157 * t29;
t46 = -t100 * t65 + t64 * t116;
t47 = t100 * t64 + t65 * t116;
t148 = t155 * t47 + t157 * t46;
t139 = m(5) + m(6) + m(7);
t130 = -t27 * t62 - t28 * t68;
t123 = t51 * t69;
t122 = t53 * t69;
t118 = t65 * t69;
t117 = t65 * t72;
t115 = t69 * t74;
t114 = t72 * t74;
t113 = -t29 * t62 - t30 * t68;
t106 = -t46 * t62 - t47 * t68;
t95 = -t27 * pkin(4) + t28 * pkin(10);
t94 = -t29 * pkin(4) + pkin(10) * t30;
t93 = -t46 * pkin(4) + pkin(10) * t47;
t88 = t151 * pkin(3);
t83 = t152 * pkin(3);
t79 = -t51 * t70 - t73 * t91;
t78 = t79 * pkin(3);
t75 = -pkin(9) - pkin(8);
t63 = pkin(3) * t73 + pkin(2);
t52 = t99 * t71 + t74 * t81;
t50 = t98 * t71 - t74 * t82;
t1 = [(-m(2) - m(3) - m(4) - t139) * g(3) (-t123 * t132 - t156 * (-t50 * t117 + t123) + t140 * (t50 * t118 + t51 * t72) - t139 * (-t50 * t63 - t51 * t75) + t138 * t51 + t133 * t50) * g(2) + (-t122 * t132 - t156 * (-t52 * t117 + t122) + t140 * (t52 * t118 + t53 * t72) - t139 * (-t52 * t63 - t53 * t75) + t138 * t53 + t133 * t52) * g(1) + ((-t139 * t63 - t133) * t74 + (-t114 * t156 - t140 * t115) * t65 + (t137 * t69 + t139 * t75 + t140 * t72 + t138) * t71) * g(3) * t67 (-t152 * mrSges(4,1) - (-t100 * t70 - t73 * t116) * mrSges(4,2) - m(5) * t83 - m(6) * (t83 + t93) - m(7) * (t83 + t106) + t148) * g(3) + (-t79 * mrSges(4,1) - (-t51 * t73 + t70 * t91) * mrSges(4,2) - m(5) * t78 - m(6) * (t78 + t95) - m(7) * (t78 + t130) + t150) * g(2) + (-t151 * mrSges(4,1) - (-t53 * t73 - t70 * t90) * mrSges(4,2) - m(5) * t88 - m(6) * (t88 + t94) - m(7) * (t88 + t113) + t149) * g(1) (-m(6) * t93 - m(7) * t106 + t148) * g(3) + (-m(6) * t95 - m(7) * t130 + t150) * g(2) + (-m(6) * t94 - m(7) * t113 + t149) * g(1) (t140 * (t67 * t115 - t47 * t72) + t137 * (-t67 * t114 - t47 * t69)) * g(3) + (t140 * (-t28 * t72 - t50 * t69) + t137 * (-t28 * t69 + t50 * t72)) * g(2) + (t140 * (-t30 * t72 - t52 * t69) + t137 * (-t30 * t69 + t52 * t72)) * g(1) (-g(1) * t29 - g(2) * t27 - g(3) * t46) * m(7)];
taug  = t1(:);
