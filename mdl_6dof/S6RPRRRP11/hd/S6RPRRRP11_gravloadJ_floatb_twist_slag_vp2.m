% Calculate Gravitation load on the joints for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:30
% EndTime: 2019-03-09 06:34:33
% DurationCPUTime: 1.50s
% Computational Cost: add. (1218->136), mult. (3285->194), div. (0->0), fcn. (4170->14), ass. (0->73)
t120 = cos(qJ(4));
t121 = cos(qJ(3));
t113 = cos(pkin(7));
t122 = cos(qJ(1));
t109 = sin(pkin(12));
t119 = sin(qJ(1));
t112 = cos(pkin(12));
t114 = cos(pkin(6));
t96 = t114 * t112;
t82 = t119 * t109 - t122 * t96;
t110 = sin(pkin(7));
t111 = sin(pkin(6));
t92 = t111 * t110;
t143 = t82 * t113 + t122 * t92;
t94 = t114 * t109;
t45 = t119 * t112 + t122 * t94;
t61 = sin(qJ(3));
t30 = -t45 * t121 + t143 * t61;
t60 = sin(qJ(4));
t93 = t113 * t111;
t69 = -t82 * t110 + t122 * t93;
t16 = t30 * t120 + t69 * t60;
t27 = t143 * t121 + t45 * t61;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t147 = t16 * t59 + t27 * t62;
t146 = t16 * t62 - t27 * t59;
t142 = -t69 * t120 + t30 * t60;
t140 = mrSges(6,1) + mrSges(7,1);
t134 = mrSges(6,2) + mrSges(7,2);
t56 = pkin(5) * t62 + pkin(4);
t138 = m(6) * pkin(4) + m(7) * t56 + mrSges(5,1);
t137 = m(6) * pkin(11) - m(7) * (-qJ(6) - pkin(11)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t75 = t122 * t109 + t119 * t96;
t63 = t75 * t110 + t119 * t93;
t123 = m(7) * pkin(5);
t135 = mrSges(4,2) - mrSges(5,3);
t133 = t75 * t113 - t119 * t92;
t132 = t114 * t110 + t112 * t93;
t130 = -m(5) - m(6) - m(7);
t129 = -t134 * t59 + t140 * t62 + t138;
t128 = -t123 - t140;
t127 = t138 * t120 + t137 * t60 + mrSges(4,1);
t125 = -m(7) * (pkin(5) * t59 + pkin(10)) + t135;
t118 = t30 * t59;
t46 = t122 * t112 - t119 * t94;
t32 = t46 * t121 - t133 * t61;
t117 = t32 * t59;
t91 = t111 * t109;
t38 = t121 * t91 + t132 * t61;
t116 = t38 * t59;
t97 = t111 * t119;
t115 = t122 * pkin(1) + qJ(2) * t97;
t106 = t59 * t120;
t105 = t62 * t120;
t98 = t122 * t111;
t100 = -t119 * pkin(1) + qJ(2) * t98;
t18 = t32 * t120 + t63 * t60;
t31 = t133 * t121 + t46 * t61;
t5 = -t18 * t59 + t31 * t62;
t74 = -t112 * t92 + t114 * t113;
t71 = -t45 * pkin(2) + t69 * pkin(9) + t100;
t68 = t30 * pkin(3) + t71;
t67 = t46 * pkin(2) + t63 * pkin(9) + t115;
t66 = t32 * pkin(3) + t67;
t65 = -pkin(10) * t27 + t68;
t64 = t31 * pkin(10) + t66;
t37 = -t132 * t121 + t61 * t91;
t26 = t38 * t120 + t74 * t60;
t25 = -t74 * t120 + t38 * t60;
t17 = -t63 * t120 + t32 * t60;
t6 = t18 * t62 + t31 * t59;
t1 = [(-t122 * mrSges(2,1) + t119 * mrSges(2,2) - m(3) * t115 - t46 * mrSges(3,1) + t75 * mrSges(3,2) - mrSges(3,3) * t97 - m(4) * t67 - t32 * mrSges(4,1) - t63 * mrSges(4,3) - m(5) * t64 - t18 * mrSges(5,1) - m(6) * (t18 * pkin(4) + t64) - m(7) * (t18 * t56 + t66) - t140 * t6 - t134 * t5 + t125 * t31 - t137 * t17) * g(2) + (t119 * mrSges(2,1) + t122 * mrSges(2,2) - m(3) * t100 + t45 * mrSges(3,1) - t82 * mrSges(3,2) - mrSges(3,3) * t98 - m(4) * t71 - t30 * mrSges(4,1) - t69 * mrSges(4,3) - m(5) * t65 - t16 * mrSges(5,1) - m(6) * (t16 * pkin(4) + t65) - m(7) * (t16 * t56 + t68) - t140 * t146 + t134 * t147 - t125 * t27 - t137 * t142) * g(1) (-g(1) * t97 + g(2) * t98 - g(3) * t114) * (m(3) + m(4) - t130) (-t116 * t123 + t135 * t38 - t140 * (-t37 * t105 + t116) - t134 * (t37 * t106 + t38 * t62) + t130 * (-t37 * pkin(3) + t38 * pkin(10)) + t127 * t37) * g(3) + (t118 * t123 - t140 * (-t27 * t105 - t118) - t134 * (t27 * t106 - t30 * t62) - t135 * t30 + t130 * (-t27 * pkin(3) - pkin(10) * t30) + t127 * t27) * g(2) + (-t117 * t123 - t134 * (t31 * t106 + t32 * t62) + t135 * t32 + t130 * (-t31 * pkin(3) + t32 * pkin(10)) - t140 * (-t31 * t105 + t117) + t127 * t31) * g(1) (t129 * t25 - t137 * t26) * g(3) + (-t129 * t142 + t137 * t16) * g(2) + (t129 * t17 - t137 * t18) * g(1) (-t134 * (-t26 * t62 - t37 * t59) + t128 * (-t26 * t59 + t37 * t62)) * g(3) + (t128 * t147 - t134 * t146) * g(2) + (t128 * t5 + t134 * t6) * g(1) (-g(1) * t17 + g(2) * t142 - g(3) * t25) * m(7)];
taug  = t1(:);
