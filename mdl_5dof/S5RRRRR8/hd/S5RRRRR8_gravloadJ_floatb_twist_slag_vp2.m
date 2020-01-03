% Calculate Gravitation load on the joints for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:23
% EndTime: 2019-12-31 22:24:25
% DurationCPUTime: 0.67s
% Computational Cost: add. (364->100), mult. (421->114), div. (0->0), fcn. (375->10), ass. (0->64)
t39 = qJ(4) + qJ(5);
t34 = sin(t39);
t41 = sin(qJ(4));
t125 = -t41 * mrSges(5,2) - t34 * mrSges(6,2);
t124 = -mrSges(5,3) - mrSges(6,3);
t36 = cos(t39);
t44 = cos(qJ(4));
t123 = t44 * mrSges(5,1) + t36 * mrSges(6,1);
t40 = qJ(2) + qJ(3);
t35 = sin(t40);
t37 = cos(t40);
t122 = t37 * (-m(5) * pkin(8) + t124) + t125 * t35;
t101 = m(6) * pkin(4);
t118 = t123 * t35;
t117 = -t37 * mrSges(4,1) + (mrSges(4,2) + t124) * t35;
t112 = t37 * pkin(3) + t35 * pkin(8);
t32 = t44 * pkin(4) + pkin(3);
t47 = -pkin(9) - pkin(8);
t114 = t37 * t32 - t35 * t47;
t116 = -m(5) * t112 - m(6) * t114;
t115 = t41 * t101;
t56 = -t32 * t35 - t37 * t47;
t97 = pkin(3) * t35;
t42 = sin(qJ(2));
t98 = pkin(2) * t42;
t111 = -m(6) * (t56 - t98) - m(5) * (-t97 - t98) + t118;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t110 = g(1) * t46 + g(2) * t43;
t109 = mrSges(5,1) + t101;
t108 = m(4) + m(5) + m(6);
t107 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t106 = t117 + (-t123 - t125) * t37;
t105 = t122 * t43;
t104 = t122 * t46;
t103 = m(5) * t97 - m(6) * t56 + t118;
t45 = cos(qJ(2));
t62 = t45 * mrSges(3,1) - t42 * mrSges(3,2);
t102 = m(3) * pkin(1) + mrSges(2,1) - t117 + t62;
t79 = t46 * t36;
t85 = t43 * t34;
t5 = t37 * t85 + t79;
t80 = t46 * t34;
t84 = t43 * t36;
t6 = -t37 * t84 + t80;
t100 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = -t37 * t80 + t84;
t8 = t37 * t79 + t85;
t99 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t93 = g(3) * t35;
t38 = t45 * pkin(2);
t83 = t43 * t41;
t82 = t43 * t44;
t78 = t46 * t41;
t77 = t46 * t44;
t59 = mrSges(4,1) * t35 + mrSges(4,2) * t37;
t58 = -mrSges(6,1) * t34 - mrSges(6,2) * t36;
t11 = -t37 * t78 + t82;
t9 = t37 * t83 + t77;
t48 = -pkin(7) - pkin(6);
t33 = t38 + pkin(1);
t12 = t37 * t77 + t83;
t10 = -t37 * t82 + t78;
t1 = [(-t83 * t101 - t12 * mrSges(5,1) - t8 * mrSges(6,1) - t11 * mrSges(5,2) - t7 * mrSges(6,2) - t108 * (t46 * t33 - t43 * t48) + t107 * t43 + (-t102 + t116) * t46) * g(2) + (-t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (t108 * t48 + t107 - t115) * t46 + (m(4) * t33 - m(5) * (-t33 - t112) - m(6) * (-t33 - t114) + t102) * t43) * g(1), (t111 * t43 + t105) * g(2) + (t111 * t46 + t104) * g(1) + (-t62 - m(4) * t38 - m(5) * (t38 + t112) - m(6) * (t38 + t114) + t106) * g(3) + (m(4) * t98 + mrSges(3,1) * t42 + mrSges(3,2) * t45 + t59) * t110, t110 * t59 + (t103 * t43 + t105) * g(2) + (t103 * t46 + t104) * g(1) + (t106 + t116) * g(3), (mrSges(5,1) * t41 + mrSges(5,2) * t44 + t115 - t58) * t93 + (-t10 * mrSges(5,2) + t109 * t9 - t100) * g(2) + (t12 * mrSges(5,2) - t109 * t11 - t99) * g(1), -g(1) * t99 - g(2) * t100 - t58 * t93];
taug = t1(:);
