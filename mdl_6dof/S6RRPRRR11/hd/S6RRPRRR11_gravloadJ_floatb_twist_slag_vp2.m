% Calculate Gravitation load on the joints for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:28:30
% EndTime: 2018-11-23 17:28:31
% DurationCPUTime: 0.90s
% Computational Cost: add. (427->138), mult. (600->154), div. (0->0), fcn. (548->10), ass. (0->65)
t111 = mrSges(3,1) - mrSges(4,2);
t110 = -mrSges(3,2) + mrSges(4,3);
t109 = -mrSges(6,3) - mrSges(7,3);
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t97 = g(1) * t50 + g(2) * t47;
t44 = qJ(4) + qJ(5);
t34 = sin(t44);
t45 = sin(qJ(4));
t87 = pkin(4) * t45;
t22 = pkin(5) * t34 + t87;
t37 = qJ(6) + t44;
t31 = sin(t37);
t32 = cos(t37);
t35 = cos(t44);
t48 = cos(qJ(4));
t108 = -m(6) * t87 - m(7) * t22 - t45 * mrSges(5,1) - t34 * mrSges(6,1) - t31 * mrSges(7,1) - t48 * mrSges(5,2) - t35 * mrSges(6,2) - t32 * mrSges(7,2);
t51 = -pkin(9) - pkin(8);
t43 = -pkin(10) + t51;
t107 = -m(5) * (-pkin(2) - pkin(8)) + mrSges(5,3) - m(6) * (-pkin(2) + t51) - m(7) * (-pkin(2) + t43) - t109;
t103 = -m(6) * pkin(4) - mrSges(5,1);
t49 = cos(qJ(2));
t24 = t49 * t31 * mrSges(7,2);
t102 = -t49 * t34 * mrSges(6,2) - t24;
t46 = sin(qJ(2));
t101 = t110 * t46 + t111 * t49;
t79 = t46 * t47;
t15 = t34 * t50 + t35 * t79;
t16 = -t34 * t79 + t35 * t50;
t7 = t31 * t50 + t32 * t79;
t8 = -t31 * t79 + t32 * t50;
t88 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t100 = -t15 * mrSges(6,1) - t16 * mrSges(6,2) - t88;
t81 = mrSges(7,1) * t32;
t99 = mrSges(6,1) * t35 + t81;
t78 = t46 * t50;
t13 = -t34 * t47 + t35 * t78;
t14 = t34 * t78 + t35 * t47;
t5 = -t31 * t47 + t32 * t78;
t6 = t31 * t78 + t32 * t47;
t89 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t98 = -t13 * mrSges(6,1) + t14 * mrSges(6,2) - t89;
t96 = m(4) + m(5) + m(6) + m(7);
t94 = t109 * t49 - t101;
t30 = pkin(5) * t35;
t39 = t48 * pkin(4);
t23 = t30 + t39;
t92 = -m(5) * pkin(3) - m(6) * (t39 + pkin(3)) - m(7) * (pkin(3) + t23) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t90 = m(7) * pkin(5);
t84 = g(3) * t49;
t40 = t49 * pkin(2);
t80 = t43 * t49;
t77 = t47 * t48;
t76 = t48 * t50;
t73 = t49 * t50;
t72 = t49 * t51;
t36 = t46 * qJ(3);
t70 = t40 + t36;
t69 = t50 * pkin(1) + t47 * pkin(7);
t66 = -pkin(1) - t36;
t17 = -t45 * t47 + t46 * t76;
t19 = t45 * t50 + t46 * t77;
t20 = -t45 * t79 + t76;
t18 = t45 * t78 + t77;
t1 = [(-m(3) * t69 - t18 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t17 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (-m(5) * pkin(8) - mrSges(5,3)) * t73 - t96 * (pkin(2) * t73 + t50 * t36 + t69) + t92 * t47 + (-mrSges(2,1) - m(6) * (t46 * t87 - t72) - m(7) * (t22 * t46 - t80) + t94) * t50) * g(2) + (-t20 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) + t19 * mrSges(5,2) + t15 * mrSges(6,2) + t7 * mrSges(7,2) + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t66 - t40) - m(5) * t66 - m(6) * (-pkin(1) + (-qJ(3) - t87) * t46) - m(7) * (-pkin(1) + (-qJ(3) - t22) * t46) + t107 * t49 + t101) * t47 + ((-m(3) - t96) * pkin(7) + t92) * t50) * g(1) (-m(4) * t70 - m(5) * (pkin(8) * t49 + t70) - t49 * mrSges(5,3) - m(6) * (t70 - t72) - m(7) * (t70 - t80) + t108 * t46 + t94) * g(3) + ((m(4) * pkin(2) + t107 + t111) * t46 + (-qJ(3) * t96 + t108 - t110) * t49) * t97 (-t97 * t46 + t84) * t96 -(-mrSges(5,1) * t48 + mrSges(5,2) * t45) * t84 + ((m(6) * t39 + m(7) * t23 + t99) * t49 + t102) * g(3) + (-t20 * mrSges(5,2) - m(7) * (t22 * t50 + t23 * t79) + t103 * t19 + t100) * g(2) + (t18 * mrSges(5,2) - m(7) * (-t22 * t47 + t23 * t78) + t103 * t17 + t98) * g(1) ((m(7) * t30 + t99) * t49 + t102) * g(3) + (-t15 * t90 + t100) * g(2) + (-t13 * t90 + t98) * g(1), -g(1) * t89 - g(2) * t88 - g(3) * (-t49 * t81 + t24)];
taug  = t1(:);
