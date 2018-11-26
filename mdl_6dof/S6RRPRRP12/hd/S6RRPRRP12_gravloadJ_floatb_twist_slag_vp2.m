% Calculate Gravitation load on the joints for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:19:05
% EndTime: 2018-11-23 17:19:06
% DurationCPUTime: 1.00s
% Computational Cost: add. (402->134), mult. (638->153), div. (0->0), fcn. (596->8), ass. (0->72)
t120 = mrSges(6,2) - mrSges(7,3);
t108 = -m(7) * qJ(6) + t120;
t121 = mrSges(6,1) + mrSges(7,1);
t109 = -m(7) * pkin(5) - t121;
t116 = -m(6) - m(7);
t119 = mrSges(7,2) + mrSges(6,3);
t45 = qJ(4) + qJ(5);
t38 = sin(t45);
t39 = cos(t45);
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t118 = -t46 * mrSges(5,1) - t49 * mrSges(5,2) - t108 * t39 + t109 * t38;
t117 = -m(4) - m(5);
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t115 = (-mrSges(3,1) + mrSges(4,2)) * t50 + (mrSges(3,2) - mrSges(4,3)) * t47;
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t114 = g(1) * t51 + g(2) * t48;
t95 = t47 * t48;
t13 = t38 * t51 + t39 * t95;
t14 = -t38 * t95 + t39 * t51;
t113 = -t120 * t14 - t121 * t13;
t94 = t47 * t51;
t11 = t38 * t48 - t39 * t94;
t12 = t38 * t94 + t39 * t48;
t112 = t121 * t11 + t120 * t12;
t69 = m(5) * (-pkin(2) - pkin(8)) - mrSges(5,3);
t111 = (-mrSges(4,3) + t118) * t50 + (-mrSges(4,2) - t69 + (m(4) - t116) * pkin(2) + t119) * t47;
t110 = -t119 * t50 + t115;
t107 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t102 = pkin(4) * t49;
t99 = g(3) * t50;
t42 = t50 * pkin(2);
t98 = t39 * t50;
t97 = t46 * t48;
t96 = t46 * t51;
t93 = t48 * t49;
t92 = t49 * t51;
t89 = t50 * t51;
t52 = -pkin(9) - pkin(8);
t88 = t50 * t52;
t87 = t51 * t52;
t35 = pkin(4) * t96;
t80 = t47 * t93;
t86 = pkin(4) * t80 + t35;
t40 = t47 * qJ(3);
t85 = t42 + t40;
t84 = t51 * pkin(1) + t48 * pkin(7);
t83 = qJ(3) * t50;
t82 = pkin(4) * t97;
t34 = t47 * t46 * pkin(4);
t81 = t46 * t94;
t79 = t47 * t92;
t74 = -pkin(1) - t40;
t73 = -t11 * pkin(5) + qJ(6) * t12;
t72 = t13 * pkin(5) - qJ(6) * t14;
t71 = pkin(2) * t89 + t51 * t40 + t84;
t70 = pkin(4) * t79 - t82;
t64 = -t39 * mrSges(7,1) - t38 * mrSges(7,3);
t62 = -pkin(5) * t39 - qJ(6) * t38;
t60 = t74 - t42;
t43 = t51 * pkin(7);
t37 = pkin(3) + t102;
t32 = t51 * t83;
t31 = t48 * t83;
t27 = t50 * t38 * mrSges(6,2);
t18 = -t46 * t95 + t92;
t17 = t80 + t96;
t16 = t81 + t93;
t15 = t79 - t97;
t1 = [(-m(3) * t84 - t16 * mrSges(5,1) - t15 * mrSges(5,2) + t117 * t71 + t116 * (pkin(4) * t81 + t48 * t37 - t50 * t87 + t71) + t109 * t12 + t108 * t11 + (-m(5) * pkin(8) - mrSges(5,3) - t119) * t89 + (-mrSges(2,1) + t115) * t51 + t107 * t48) * g(2) + (-t18 * mrSges(5,1) + t17 * mrSges(5,2) + t116 * (t51 * t37 + t48 * t88 + t43) + t109 * t14 + t108 * t13 + (-m(3) + t117) * t43 + t107 * t51 + (m(3) * pkin(1) - m(4) * t60 - m(5) * t74 - t69 * t50 + mrSges(2,1) + t116 * (t60 - t34) - t110) * t48) * g(1), t114 * (mrSges(3,1) * t47 + mrSges(3,2) * t50) + (t116 * (t50 * t82 + t52 * t95 + t31) + t117 * t31 + t111 * t48) * g(2) + (t116 * (t50 * t35 + t47 * t87 + t32) + t117 * t32 + t111 * t51) * g(1) + (-m(4) * t85 - m(5) * (pkin(8) * t50 + t85) - t50 * mrSges(5,3) + t116 * (t34 + t85 - t88) + t118 * t47 + t110) * g(3) (-t114 * t47 + t99) * (-t116 - t117) -g(3) * (t27 + (-m(6) * t102 - mrSges(6,1) * t39) * t50) + (mrSges(5,1) * t49 - mrSges(5,2) * t46 - m(7) * (t62 - t102) - t64) * t99 + (-t17 * mrSges(5,1) - t18 * mrSges(5,2) - m(6) * t86 - m(7) * (t72 + t86) + t113) * g(2) + (-t15 * mrSges(5,1) + t16 * mrSges(5,2) - m(6) * t70 - m(7) * (t70 + t73) + t112) * g(1), -g(3) * (-mrSges(6,1) * t98 + t27) - (m(7) * t62 + t64) * t99 + (-m(7) * t72 + t113) * g(2) + (-m(7) * t73 + t112) * g(1) (-g(1) * t11 + g(2) * t13 - g(3) * t98) * m(7)];
taug  = t1(:);
