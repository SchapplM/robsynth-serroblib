% Calculate Gravitation load on the joints for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:11:16
% EndTime: 2018-11-23 17:11:17
% DurationCPUTime: 0.82s
% Computational Cost: add. (551->109), mult. (527->122), div. (0->0), fcn. (469->10), ass. (0->62)
t118 = mrSges(6,2) - mrSges(7,3);
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t117 = -t46 * mrSges(7,1) + t118 * t43;
t116 = mrSges(6,3) + mrSges(7,2);
t41 = qJ(2) + pkin(10);
t38 = qJ(4) + t41;
t31 = sin(t38);
t77 = qJ(6) * t43;
t114 = (mrSges(6,1) * t46 - m(7) * (-pkin(5) * t46 - pkin(4) - t77) - t117) * t31;
t107 = m(6) + m(7);
t113 = m(5) + t107;
t112 = t116 * t31;
t36 = sin(t41);
t37 = cos(t41);
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t111 = -t47 * mrSges(3,1) - t37 * mrSges(4,1) + t44 * mrSges(3,2) + t36 * mrSges(4,2);
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t108 = g(1) * t48 + g(2) * t45;
t32 = cos(t38);
t89 = t32 * t46;
t106 = pkin(5) * t89 + t32 * t77;
t105 = t32 * mrSges(5,1) - t31 * mrSges(5,2);
t79 = t32 * pkin(4) + t31 * pkin(9);
t88 = t32 * t48;
t104 = t114 * t48 - t116 * t88;
t90 = t32 * t45;
t103 = t114 * t45 - t116 * t90;
t102 = -mrSges(6,1) * t89 + t117 * t32 - t105 - t112;
t39 = t47 * pkin(2);
t100 = -mrSges(2,1) - m(3) * pkin(1) - m(4) * (t39 + pkin(1)) - t105 + t111;
t42 = -qJ(3) - pkin(7);
t99 = -m(3) * pkin(7) + m(4) * t42 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t98 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t97 = m(7) * qJ(6) - t118;
t96 = pkin(2) * t44;
t93 = g(3) * t31;
t91 = t31 * t48;
t84 = t45 * t43;
t83 = t45 * t46;
t81 = t46 * t48;
t80 = t48 * t43;
t78 = pkin(3) * t37 + t39;
t10 = pkin(1) + t78;
t40 = -pkin(8) + t42;
t70 = t48 * t10 - t40 * t45;
t68 = t78 + t79;
t20 = pkin(9) * t90;
t67 = -pkin(4) * t31 * t45 + t20;
t23 = pkin(9) * t88;
t66 = -pkin(4) * t91 + t23;
t61 = mrSges(5,1) * t31 + mrSges(5,2) * t32;
t11 = -pkin(3) * t36 - t96;
t7 = t48 * t11;
t6 = t45 * t11;
t4 = t32 * t81 + t84;
t3 = t32 * t80 - t83;
t2 = t32 * t83 - t80;
t1 = t32 * t84 + t81;
t5 = [(-m(5) * t70 - t116 * t91 - t107 * (pkin(4) * t88 + pkin(9) * t91 + t70) - t98 * t4 - t97 * t3 + t100 * t48 + t99 * t45) * g(2) + (t98 * t2 + t97 * t1 + (m(5) * t10 - t107 * (-t10 - t79) - t100 + t112) * t45 + (t113 * t40 + t99) * t48) * g(1) (-m(6) * (t6 + t67) - m(7) * (t20 + t6) + t103) * g(2) + (-m(6) * (t66 + t7) - m(7) * (t23 + t7) + t104) * g(1) + (-m(4) * t39 - m(5) * t78 - m(6) * t68 - m(7) * (t68 + t106) + t102 + t111) * g(3) + t108 * (m(4) * t96 - m(5) * t11 + mrSges(3,1) * t44 + mrSges(4,1) * t36 + mrSges(3,2) * t47 + mrSges(4,2) * t37 + t61) (-g(1) * t45 + g(2) * t48) * (m(4) + t113) t108 * t61 + (-m(6) * t67 - m(7) * t20 + t103) * g(2) + (-m(6) * t66 - m(7) * t23 + t104) * g(1) + (-m(6) * t79 - m(7) * (t79 + t106) + t102) * g(3) (t98 * t43 - t97 * t46) * t93 + (t98 * t1 - t97 * t2) * g(2) + (t98 * t3 - t97 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t43 * t93) * m(7)];
taug  = t5(:);
