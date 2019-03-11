% Calculate Gravitation load on the joints for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:12:02
% EndTime: 2019-03-09 22:12:05
% DurationCPUTime: 1.23s
% Computational Cost: add. (638->150), mult. (866->175), div. (0->0), fcn. (864->10), ass. (0->75)
t150 = mrSges(6,2) + mrSges(5,3);
t141 = m(6) + m(7);
t49 = qJ(2) + qJ(3);
t46 = sin(t49);
t51 = sin(qJ(4));
t117 = t46 * t51;
t50 = sin(qJ(6));
t54 = cos(qJ(6));
t55 = cos(qJ(4));
t76 = t50 * t51 + t54 * t55;
t140 = t76 * t46;
t77 = t50 * t55 - t51 * t54;
t9 = t77 * t46;
t148 = mrSges(7,1) * t140 - mrSges(5,2) * t117 - mrSges(7,2) * t9;
t53 = sin(qJ(1));
t111 = t53 * t51;
t47 = cos(t49);
t57 = cos(qJ(1));
t16 = t111 * t47 + t55 * t57;
t107 = t57 * t51;
t110 = t53 * t55;
t17 = t110 * t47 - t107;
t135 = t16 * t54 - t17 * t50;
t78 = t16 * t50 + t17 * t54;
t147 = t135 * mrSges(7,1) - t78 * mrSges(7,2);
t146 = t150 * t46;
t144 = t55 * mrSges(6,1) + t51 * mrSges(6,3);
t143 = g(1) * t57 + g(2) * t53;
t116 = t46 * t55;
t124 = -pkin(4) - pkin(5);
t102 = qJ(5) * t51;
t93 = -pkin(3) - t102;
t142 = mrSges(5,1) * t116 + (-m(6) * (-pkin(4) * t55 + t93) + t144 - m(7) * (t124 * t55 + t93)) * t46;
t138 = mrSges(5,1) + mrSges(6,1);
t137 = mrSges(5,2) - mrSges(6,3);
t134 = t47 * mrSges(4,1) - t46 * mrSges(4,2);
t52 = sin(qJ(2));
t56 = cos(qJ(2));
t82 = t56 * mrSges(3,1) - t52 * mrSges(3,2);
t133 = -m(3) * pkin(1) - mrSges(2,1) - t134 - t82;
t132 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t122 = pkin(3) * t46;
t123 = pkin(2) * t52;
t131 = -m(7) * (-pkin(10) * t47 - t123) + t47 * mrSges(7,3) + m(6) * t123 - m(5) * (-t122 - t123) + t142;
t128 = -m(7) * pkin(10) - mrSges(7,3);
t130 = -t128 * t47 + t142;
t108 = t57 * t47;
t36 = pkin(9) * t108;
t129 = -t108 * t150 - t141 * t36 + t148 * t57;
t114 = t47 * t55;
t127 = -mrSges(5,1) * t114 + t46 * mrSges(7,3) - t134 - t146 + (-mrSges(7,1) * t76 + t51 * mrSges(5,2) + mrSges(7,2) * t77 - t144) * t47;
t126 = m(7) * pkin(5) + t138;
t125 = (t148 + (-t150 + (-m(5) - t141) * pkin(9)) * t47) * t53;
t41 = t46 * pkin(9);
t42 = t47 * pkin(3);
t48 = t56 * pkin(2);
t115 = t46 * t57;
t58 = -pkin(8) - pkin(7);
t106 = t57 * t58;
t103 = t42 + t41;
t45 = t48 + pkin(1);
t95 = -t45 - t42;
t91 = t57 * t45 - t53 * t58;
t88 = pkin(4) * t114 + t47 * t102 + t103;
t18 = t107 * t47 - t110;
t19 = t108 * t55 + t111;
t1 = t18 * t54 - t19 * t50;
t2 = t18 * t50 + t19 * t54;
t84 = mrSges(7,1) * t1 - mrSges(7,2) * t2;
t83 = -t9 * mrSges(7,1) - mrSges(7,2) * t140;
t79 = mrSges(4,1) * t46 + mrSges(4,2) * t47;
t74 = pkin(3) * t108 + pkin(9) * t115 + t91;
t71 = pkin(5) * t114 - pkin(10) * t46 + t88;
t22 = qJ(5) * t116;
t3 = [(-m(4) * t91 - m(5) * t74 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t141 * (t19 * pkin(4) + t18 * qJ(5) + t74) - t126 * t19 + t137 * t18 + t133 * t57 + t132 * t53 + (-t128 - t150) * t115) * g(2) + (m(5) * t106 + t78 * mrSges(7,1) + t135 * mrSges(7,2) - t141 * (-t17 * pkin(4) - t16 * qJ(5) - t106) + t126 * t17 - t137 * t16 + (m(4) * t58 + t132) * t57 + (m(4) * t45 - m(7) * t95 - (m(7) * (-pkin(9) + pkin(10)) + mrSges(7,3)) * t46 + (-m(5) - m(6)) * (t95 - t41) - t133 + t146) * t53) * g(1) (t131 * t53 + t125) * g(2) + (-m(5) * t36 + t131 * t57 + t129) * g(1) + (-t82 - m(4) * t48 - m(5) * (t48 + t103) - m(6) * (t48 + t88) - m(7) * (t48 + t71) + t127) * g(3) + t143 * (m(4) * t123 + mrSges(3,1) * t52 + mrSges(3,2) * t56 + t79) t143 * t79 + ((m(5) * t122 + t130) * t53 + t125) * g(2) + (-m(5) * (-pkin(3) * t115 + t36) + t130 * t57 + t129) * g(1) + (-m(5) * t103 - m(6) * t88 - m(7) * t71 + t127) * g(3) (-m(6) * t22 - m(7) * (t117 * t124 + t22) + t83 + (t137 * t55 + (m(6) * pkin(4) + t138) * t51) * t46) * g(3) + (-t141 * (-t16 * pkin(4) + qJ(5) * t17) + t137 * t17 + t126 * t16 + t147) * g(2) + (t84 - t141 * (-t18 * pkin(4) + qJ(5) * t19) + t137 * t19 + t126 * t18) * g(1), t141 * (-g(1) * t18 - g(2) * t16 - g(3) * t117) -g(1) * t84 - g(2) * t147 - g(3) * t83];
taug  = t3(:);
