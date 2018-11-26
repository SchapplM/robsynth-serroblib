% Calculate Gravitation load on the joints for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2018-11-23 16:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:50:54
% EndTime: 2018-11-23 16:50:55
% DurationCPUTime: 1.00s
% Computational Cost: add. (1377->128), mult. (1159->166), div. (0->0), fcn. (1056->20), ass. (0->76)
t109 = m(6) + m(7);
t118 = m(5) + t109;
t114 = qJ(4) * t118 - mrSges(4,2) + mrSges(5,3);
t115 = -t109 * pkin(9) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t60 = sin(qJ(6));
t64 = cos(qJ(6));
t113 = t60 * mrSges(7,1) + t64 * mrSges(7,2) - t115;
t117 = m(7) * pkin(5) + t64 * mrSges(7,1) - t60 * mrSges(7,2) + mrSges(6,1);
t119 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t56 = pkin(6) + qJ(2);
t47 = cos(t56) / 0.2e1;
t57 = pkin(6) - qJ(2);
t54 = cos(t57);
t36 = t54 / 0.2e1 + t47;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t112 = -t117 * t61 + t119 * t65 - t114;
t111 = sin(t56) / 0.2e1;
t51 = sin(t57);
t108 = pkin(2) * t51;
t58 = sin(pkin(6));
t63 = sin(qJ(1));
t107 = t58 * t63;
t67 = cos(qJ(1));
t106 = t58 * t67;
t62 = sin(qJ(2));
t105 = t62 * t67;
t55 = qJ(2) + pkin(11);
t52 = cos(t55);
t104 = t63 * t52;
t103 = t63 * t62;
t102 = t67 * t52;
t90 = pkin(6) + t55;
t83 = sin(t90);
t78 = t83 / 0.2e1;
t91 = pkin(6) - t55;
t84 = sin(t91);
t72 = t78 - t84 / 0.2e1;
t20 = -t63 * t72 + t102;
t66 = cos(qJ(2));
t48 = pkin(2) * t66 + pkin(1);
t37 = t67 * t48;
t99 = t20 * pkin(3) + t37;
t41 = pkin(2) * t111;
t98 = t108 / 0.2e1 + t41;
t95 = pkin(4) * t107 + t99;
t93 = m(4) + t118;
t34 = t36 * pkin(2);
t88 = -pkin(2) * t103 + t67 * t34;
t85 = cos(t90);
t82 = m(3) * pkin(1) + t66 * mrSges(3,1) + mrSges(2,1);
t81 = -pkin(2) * t105 - t63 * t34;
t49 = sin(t55);
t80 = cos(t91) / 0.2e1;
t70 = t85 / 0.2e1 + t80;
t16 = t49 * t63 - t67 * t70;
t7 = t65 * t106 - t16 * t61;
t5 = t61 * t106 + t16 * t65;
t79 = t84 / 0.2e1;
t73 = -t83 / 0.2e1 + t79;
t35 = t111 - t51 / 0.2e1;
t69 = t35 * mrSges(3,1) + mrSges(2,2) + (-m(3) * pkin(8) - mrSges(5,1) - mrSges(3,3) - mrSges(4,3)) * t58 + t93 * (-t108 / 0.2e1 + t41 - t58 * (pkin(8) + qJ(3)));
t59 = cos(pkin(6));
t33 = t80 - t85 / 0.2e1;
t32 = t79 + t78;
t25 = -t63 * t36 - t105;
t24 = -t36 * t67 + t103;
t23 = -t32 * t61 + t59 * t65;
t19 = t67 * t49 + t63 * t70;
t17 = t67 * t72 + t104;
t13 = t17 * pkin(3);
t4 = t65 * t107 + t19 * t61;
t3 = t61 * t107 - t19 * t65;
t2 = t20 * t60 + t4 * t64;
t1 = t20 * t64 - t4 * t60;
t6 = [(-t25 * mrSges(3,2) - m(4) * t37 - m(5) * t99 - m(6) * t95 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t95) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t119 * t3 - t114 * t19 - t82 * t67 + t115 * t20 + t69 * t63) * g(2) + (-t24 * mrSges(3,2) + m(5) * t13 - t119 * t5 + t114 * t16 - t117 * t7 + t113 * t17 + (t93 * t48 + t82) * t63 + t69 * t67 + t109 * (-pkin(4) * t106 + t13)) * g(1) (-(t111 + t51 / 0.2e1) * mrSges(3,1) - (t47 - t54 / 0.2e1) * mrSges(3,2) - m(4) * t98 - t118 * (t32 * pkin(3) + t98) + t112 * t33 - t113 * t32) * g(3) + (t24 * mrSges(3,1) - (-t35 * t67 - t63 * t66) * mrSges(3,2) - m(4) * t88 - t118 * (-t16 * pkin(3) + t88) + t112 * (-t67 * t73 + t104) + t113 * t16) * g(2) + (-t25 * mrSges(3,1) - (t63 * t35 - t66 * t67) * mrSges(3,2) - m(4) * t81 - t118 * (-t19 * pkin(3) + t81) + t112 * (t63 * t73 + t102) + t113 * t19) * g(1) (-t59 * g(3) + (-t63 * g(1) + t67 * g(2)) * t58) * t93, t118 * (-g(1) * t19 - g(2) * t16 + g(3) * t32) (-t119 * t23 - t117 * (-t32 * t65 - t59 * t61)) * g(3) + (-t117 * t5 + t119 * t7) * g(2) + (t117 * t3 - t119 * t4) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t17 * t64 + t60 * t7) * mrSges(7,1) + (-t17 * t60 + t64 * t7) * mrSges(7,2)) - g(3) * ((-t23 * t60 + t33 * t64) * mrSges(7,1) + (-t23 * t64 - t33 * t60) * mrSges(7,2))];
taug  = t6(:);
