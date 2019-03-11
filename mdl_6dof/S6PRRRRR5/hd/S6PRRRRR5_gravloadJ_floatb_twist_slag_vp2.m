% Calculate Gravitation load on the joints for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:04:01
% EndTime: 2019-03-09 01:04:05
% DurationCPUTime: 1.40s
% Computational Cost: add. (1275->136), mult. (3352->207), div. (0->0), fcn. (4249->16), ass. (0->75)
t138 = -m(5) - m(6);
t59 = qJ(5) + qJ(6);
t57 = sin(t59);
t58 = cos(t59);
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t140 = mrSges(5,1) + m(7) * (pkin(5) * t63 + pkin(4)) + t58 * mrSges(7,1) - t57 * mrSges(7,2) + m(6) * pkin(4) + t63 * mrSges(6,1) - t60 * mrSges(6,2);
t126 = mrSges(5,2) + m(7) * (-pkin(12) - pkin(11)) - mrSges(7,3) - m(6) * pkin(11) - mrSges(6,3);
t131 = m(7) - t138;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t139 = pkin(3) * t131 - t126 * t61 + t140 * t64 + mrSges(4,1);
t125 = -t60 * mrSges(6,1) - t63 * mrSges(6,2) - m(7) * (pkin(5) * t60 + pkin(10)) - t57 * mrSges(7,1) - t58 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t113 = cos(pkin(13));
t114 = cos(pkin(7));
t110 = sin(pkin(13));
t117 = sin(qJ(2));
t119 = cos(qJ(2));
t115 = cos(pkin(6));
t92 = t115 * t113;
t70 = t110 * t117 - t119 * t92;
t111 = sin(pkin(7));
t112 = sin(pkin(6));
t88 = t112 * t111;
t136 = t113 * t88 + t70 * t114;
t90 = t115 * t110;
t71 = t113 * t117 + t119 * t90;
t87 = t112 * t110;
t135 = -t111 * t87 + t71 * t114;
t89 = t114 * t112;
t134 = t115 * t111 + t119 * t89;
t133 = -m(7) * pkin(5) - mrSges(6,1);
t132 = -t111 * mrSges(4,3) + mrSges(3,2);
t127 = pkin(10) * t138 + t125;
t118 = cos(qJ(3));
t46 = t110 * t119 + t117 * t92;
t62 = sin(qJ(3));
t17 = t136 * t118 + t46 * t62;
t18 = t46 * t118 - t136 * t62;
t34 = t111 * t70 - t113 * t89;
t8 = t18 * t64 + t34 * t61;
t122 = (t17 * t58 - t57 * t8) * mrSges(7,1) + (-t17 * t57 - t58 * t8) * mrSges(7,2);
t47 = t113 * t119 - t117 * t90;
t20 = t47 * t118 - t135 * t62;
t35 = t111 * t71 + t114 * t87;
t10 = t20 * t64 + t35 * t61;
t19 = t135 * t118 + t47 * t62;
t121 = (-t10 * t57 + t19 * t58) * mrSges(7,1) + (-t10 * t58 - t19 * t57) * mrSges(7,2);
t97 = t112 * t117;
t33 = t118 * t97 + t134 * t62;
t45 = t114 * t115 - t119 * t88;
t22 = t33 * t64 + t45 * t61;
t32 = -t134 * t118 + t62 * t97;
t120 = (-t22 * t57 + t32 * t58) * mrSges(7,1) + (-t22 * t58 - t32 * t57) * mrSges(7,2);
t77 = t117 * t88;
t98 = t119 * t112;
t116 = pkin(2) * t98 + pkin(9) * t77;
t79 = t117 * t89;
t42 = t118 * t98 - t62 * t79;
t109 = t42 * pkin(3) + t116;
t105 = pkin(9) * t111;
t104 = t61 * t111;
t103 = t62 * t114;
t102 = t64 * t111;
t99 = t114 * t118;
t86 = -t70 * pkin(2) + t105 * t46;
t85 = -t71 * pkin(2) + t105 * t47;
t26 = -t103 * t46 - t118 * t70;
t84 = t26 * pkin(3) + t86;
t28 = -t103 * t47 - t118 * t71;
t83 = t28 * pkin(3) + t85;
t41 = t118 * t79 + t62 * t98;
t27 = t47 * t99 - t62 * t71;
t25 = t46 * t99 - t62 * t70;
t1 = [(-m(2) - m(3) - m(4) - t131) * g(3) (-m(4) * t116 - m(7) * t109 - mrSges(3,1) * t98 - t42 * mrSges(4,1) + mrSges(3,2) * t97 - mrSges(4,3) * t77 + t138 * (t41 * pkin(10) + t109) - t140 * (t42 * t64 + t61 * t77) + t125 * t41 + t126 * (t42 * t61 - t64 * t77)) * g(3) + (-m(4) * t86 - m(7) * t84 + t70 * mrSges(3,1) - t26 * mrSges(4,1) + t138 * (t25 * pkin(10) + t84) + t132 * t46 - t140 * (t104 * t46 + t26 * t64) + t125 * t25 + t126 * (-t102 * t46 + t26 * t61)) * g(2) + (-m(4) * t85 - m(7) * t83 + t71 * mrSges(3,1) - t28 * mrSges(4,1) + t138 * (t27 * pkin(10) + t83) + t132 * t47 - t140 * (t104 * t47 + t28 * t64) + t125 * t27 + t126 * (-t102 * t47 + t28 * t61)) * g(1) (t127 * t33 + t139 * t32) * g(3) + (t127 * t18 + t139 * t17) * g(2) + (t127 * t20 + t139 * t19) * g(1) (t126 * t22 - t140 * (-t33 * t61 + t45 * t64)) * g(3) + (t126 * t8 - t140 * (-t18 * t61 + t34 * t64)) * g(2) + (-t140 * (-t20 * t61 + t35 * t64) + t126 * t10) * g(1) (-(-t22 * t63 - t32 * t60) * mrSges(6,2) - t120 + t133 * (-t22 * t60 + t32 * t63)) * g(3) + (-(-t17 * t60 - t63 * t8) * mrSges(6,2) - t122 + t133 * (t17 * t63 - t60 * t8)) * g(2) + (-(-t10 * t63 - t19 * t60) * mrSges(6,2) - t121 + t133 * (-t10 * t60 + t19 * t63)) * g(1), -g(1) * t121 - g(2) * t122 - g(3) * t120];
taug  = t1(:);
