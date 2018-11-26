% Calculate Gravitation load on the joints for
% S6PRRRRP3
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:28:52
% EndTime: 2018-11-23 15:28:53
% DurationCPUTime: 1.06s
% Computational Cost: add. (1363->121), mult. (1537->163), div. (0->0), fcn. (1517->16), ass. (0->69)
t59 = qJ(4) + qJ(5);
t56 = cos(t59);
t63 = cos(qJ(4));
t57 = t63 * pkin(4);
t46 = pkin(5) * t56 + t57;
t60 = sin(qJ(4));
t132 = -m(6) * (t57 + pkin(3)) - m(7) * (pkin(3) + t46) - mrSges(4,1) - m(5) * pkin(3) - t63 * mrSges(5,1) + t60 * mrSges(5,2);
t66 = -pkin(10) - pkin(9);
t120 = m(6) * t66 + m(7) * (-qJ(6) + t66) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(9) - mrSges(5,3);
t131 = mrSges(6,1) + mrSges(7,1);
t127 = -mrSges(6,2) - mrSges(7,2);
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t130 = t120 * t61 + t132 * t64 - mrSges(3,1);
t126 = -m(6) * pkin(4) - mrSges(5,1);
t104 = cos(pkin(6));
t99 = pkin(6) + qJ(2);
t87 = cos(t99) / 0.2e1;
t100 = pkin(6) - qJ(2);
t90 = cos(t100);
t44 = t87 - t90 / 0.2e1;
t38 = t104 * t61 - t44 * t64;
t88 = sin(t99);
t85 = t88 / 0.2e1;
t89 = sin(t100);
t86 = t89 / 0.2e1;
t43 = t85 + t86;
t55 = sin(t59);
t21 = -t38 * t55 - t43 * t56;
t125 = t127 * (-t38 * t56 + t43 * t55) - t131 * t21;
t101 = sin(pkin(11));
t75 = t85 - t89 / 0.2e1;
t103 = cos(pkin(11));
t65 = cos(qJ(2));
t92 = t103 * t65;
t35 = -t101 * t75 + t92;
t102 = sin(pkin(6));
t77 = t102 * t101;
t28 = t35 * t64 + t61 * t77;
t62 = sin(qJ(2));
t68 = t90 / 0.2e1 + t87;
t34 = t101 * t68 + t103 * t62;
t11 = -t28 * t55 + t34 * t56;
t124 = t127 * (-t28 * t56 - t34 * t55) - t131 * t11;
t91 = t101 * t65;
t32 = t103 * t75 + t91;
t78 = t103 * t102;
t26 = t32 * t64 - t61 * t78;
t31 = t101 * t62 - t103 * t68;
t9 = -t26 * t55 + t31 * t56;
t123 = -t131 * t9 + t127 * (-t26 * t56 - t31 * t55);
t122 = -m(4) - m(6) - m(7);
t121 = t127 * t55 + t131 * t56 - t132;
t111 = pkin(4) * t60;
t45 = pkin(5) * t55 + t111;
t118 = m(5) * pkin(8) + m(6) * t111 + m(7) * t45 + t60 * mrSges(5,1) + t63 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3);
t116 = m(7) * pkin(5);
t110 = t55 * t64;
t109 = t56 * t64;
t76 = t86 - t88 / 0.2e1;
t41 = t43 * pkin(2);
t37 = -t104 * t64 - t44 * t61;
t36 = t101 * t76 + t92;
t33 = -t103 * t76 + t91;
t30 = t34 * pkin(2);
t29 = t31 * pkin(2);
t27 = t35 * t61 - t64 * t77;
t25 = t32 * t61 + t64 * t78;
t1 = [(-m(2) - m(3) - m(5) + t122) * g(3) (-m(5) * t41 - t131 * (t43 * t109 - t44 * t55) + t127 * (-t43 * t110 - t44 * t56) + t122 * (-t44 * pkin(8) + t41) + t118 * t44 + t130 * t43) * g(3) + (m(5) * t29 - t131 * (-t31 * t109 + t33 * t55) + t127 * (t31 * t110 + t33 * t56) + t122 * (t33 * pkin(8) - t29) - t118 * t33 - t130 * t31) * g(2) + (m(5) * t30 - t131 * (-t34 * t109 + t36 * t55) + t127 * (t34 * t110 + t36 * t56) + t122 * (t36 * pkin(8) - t30) - t118 * t36 - t130 * t34) * g(1) (t120 * t38 + t121 * t37) * g(3) + (t120 * t26 + t121 * t25) * g(2) + (t120 * t28 + t121 * t27) * g(1) (-(-t38 * t63 + t43 * t60) * mrSges(5,2) - m(7) * (-t38 * t45 - t43 * t46) + t126 * (-t38 * t60 - t43 * t63) + t125) * g(3) + (-(-t26 * t63 - t31 * t60) * mrSges(5,2) - m(7) * (-t26 * t45 + t31 * t46) + t126 * (-t26 * t60 + t31 * t63) + t123) * g(2) + (-(-t28 * t63 - t34 * t60) * mrSges(5,2) - m(7) * (-t28 * t45 + t34 * t46) + t126 * (-t28 * t60 + t34 * t63) + t124) * g(1) (-t21 * t116 + t125) * g(3) + (-t9 * t116 + t123) * g(2) + (-t11 * t116 + t124) * g(1) (-g(1) * t27 - g(2) * t25 - g(3) * t37) * m(7)];
taug  = t1(:);
