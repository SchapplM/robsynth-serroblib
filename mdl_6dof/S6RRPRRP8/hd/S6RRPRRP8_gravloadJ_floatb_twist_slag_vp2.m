% Calculate Gravitation load on the joints for
% S6RRPRRP8
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:00
% EndTime: 2019-03-09 12:22:02
% DurationCPUTime: 0.96s
% Computational Cost: add. (579->129), mult. (635->145), div. (0->0), fcn. (597->10), ass. (0->67)
t116 = mrSges(6,1) + mrSges(7,1);
t115 = mrSges(6,2) - mrSges(7,3);
t117 = m(6) + m(7);
t114 = -mrSges(6,3) - mrSges(7,2);
t42 = pkin(10) + qJ(4);
t36 = qJ(5) + t42;
t31 = sin(t36);
t32 = cos(t36);
t44 = cos(pkin(10));
t33 = t44 * pkin(3) + pkin(2);
t34 = sin(t42);
t35 = cos(t42);
t113 = m(5) * t33 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t115 * t31 + t116 * t32;
t46 = sin(qJ(2));
t112 = t46 * mrSges(5,3) + mrSges(2,1);
t47 = sin(qJ(1));
t48 = cos(qJ(2));
t49 = cos(qJ(1));
t79 = t48 * t49;
t17 = -t34 * t79 + t47 * t35;
t111 = g(1) * t49 + g(2) * t47;
t110 = -m(7) * qJ(6) - mrSges(7,3);
t109 = -m(4) - m(5);
t43 = sin(pkin(10));
t58 = m(4) * pkin(2) + t44 * mrSges(4,1) - t43 * mrSges(4,2);
t106 = t58 * t48;
t78 = t49 * t31;
t13 = -t47 * t32 + t48 * t78;
t14 = t47 * t31 + t32 * t79;
t105 = t115 * t14 + t116 * t13;
t80 = t47 * t48;
t11 = t31 * t80 + t32 * t49;
t12 = t32 * t80 - t78;
t104 = t116 * t11 + t115 * t12;
t103 = -m(3) + t109;
t66 = t48 * mrSges(3,1) - t46 * mrSges(3,2);
t102 = t114 * t46 - t66;
t101 = m(7) * pkin(5) + t116;
t100 = -mrSges(6,2) - t110;
t93 = pkin(3) * t43;
t98 = -m(5) * t93 - mrSges(4,1) * t43 - t44 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t92 = pkin(4) * t34;
t89 = g(3) * t46;
t88 = mrSges(6,2) * t32;
t45 = -pkin(8) - qJ(3);
t41 = -pkin(9) + t45;
t86 = t41 * t46;
t82 = t46 * t49;
t23 = pkin(4) * t35 + t33;
t20 = t48 * t23;
t77 = t49 * pkin(1) + t47 * pkin(7);
t73 = -m(5) * t45 + mrSges(5,3);
t72 = t110 * t32 * t46;
t70 = m(4) * qJ(3) + mrSges(4,3);
t69 = -t11 * pkin(5) + t12 * qJ(6);
t67 = -t13 * pkin(5) + t14 * qJ(6);
t61 = pkin(5) * t32 + qJ(6) * t31;
t60 = t33 * t48 - t45 * t46;
t59 = t17 * pkin(4);
t15 = t34 * t80 + t35 * t49;
t56 = t15 * pkin(4);
t54 = t70 * t46 + t106;
t39 = t49 * pkin(7);
t24 = t92 + t93;
t18 = t47 * t34 + t35 * t79;
t16 = t34 * t49 - t35 * t80;
t1 = [(-t18 * mrSges(5,1) - t17 * mrSges(5,2) + t114 * t82 - t117 * (t23 * t79 + t47 * t24 - t41 * t82 + t77) - t101 * t14 - t100 * t13 + t103 * t77 + t98 * t47 + (-m(5) * t60 - t112 - t54 - t66) * t49) * g(2) + (-t16 * mrSges(5,1) - t15 * mrSges(5,2) - t117 * (t49 * t24 + t47 * t86 + t39) + t101 * t12 + t100 * t11 + t103 * t39 + t98 * t49 + (m(3) * pkin(1) - m(4) * (-qJ(3) * t46 - pkin(1)) + t46 * mrSges(4,3) + t106 - m(5) * (-pkin(1) - t60) - t117 * (-pkin(1) - t20) - t102 + t112) * t47) * g(1) (-t54 - t117 * (t20 - t86) + t102) * g(3) + ((-m(7) * t61 - t113) * g(3) + t111 * (t117 * t41 + mrSges(3,2) + t114 - t70 - t73)) * t48 + (-t73 * g(3) + t111 * (mrSges(3,1) + t58 + m(6) * t23 - m(7) * (-t23 - t61) + t113)) * t46 (t48 * g(3) - t111 * t46) * (t117 - t109) -g(3) * ((m(7) * (-pkin(5) * t31 - t92) - t31 * mrSges(7,1)) * t46 - t72) + (m(6) * t92 + mrSges(5,1) * t34 + mrSges(6,1) * t31 + mrSges(5,2) * t35 + t88) * t89 + (t15 * mrSges(5,1) - t16 * mrSges(5,2) + m(6) * t56 - m(7) * (-t56 + t69) + t104) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t59 - m(7) * (t59 + t67) + t105) * g(1) ((t101 * t31 + t88) * t46 + t72) * g(3) + (-m(7) * t69 + t104) * g(2) + (-m(7) * t67 + t105) * g(1) (-g(1) * t13 - g(2) * t11 - t31 * t89) * m(7)];
taug  = t1(:);
