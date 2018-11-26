% Calculate Gravitation load on the joints for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2018-11-23 14:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:49:17
% EndTime: 2018-11-23 14:49:18
% DurationCPUTime: 0.95s
% Computational Cost: add. (2192->117), mult. (1914->167), div. (0->0), fcn. (1814->28), ass. (0->82)
t125 = -m(6) - m(7);
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t126 = m(7) * pkin(5) + t70 * mrSges(7,1) - t67 * mrSges(7,2) + mrSges(6,1);
t122 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t62 = sin(pkin(6));
t64 = cos(pkin(11));
t113 = t62 * t64;
t58 = qJ(3) + pkin(13);
t105 = pkin(7) + t58;
t89 = sin(t105) / 0.2e1;
t106 = pkin(7) - t58;
t96 = sin(t106);
t123 = t89 - t96 / 0.2e1;
t112 = sin(pkin(11));
t110 = pkin(6) - pkin(12);
t102 = sin(t110);
t109 = pkin(6) + pkin(12);
t94 = sin(t109) / 0.2e1;
t42 = t94 - t102 / 0.2e1;
t63 = cos(pkin(12));
t30 = t112 * t63 + t64 * t42;
t90 = cos(t106) / 0.2e1;
t97 = cos(t105);
t39 = t90 - t97 / 0.2e1;
t55 = cos(t58);
t111 = sin(pkin(12));
t103 = cos(t109);
t95 = cos(t110) / 0.2e1;
t81 = t95 + t103 / 0.2e1;
t73 = t112 * t111 - t64 * t81;
t10 = -t39 * t113 - t123 * t73 + t30 * t55;
t107 = t62 * t112;
t31 = -t112 * t42 + t64 * t63;
t74 = t64 * t111 + t112 * t81;
t13 = t39 * t107 - t123 * t74 + t31 * t55;
t43 = t95 - t103 / 0.2e1;
t66 = cos(pkin(6));
t80 = t94 + t102 / 0.2e1;
t17 = t123 * t80 + t66 * t39 + t43 * t55;
t59 = pkin(7) + qJ(3);
t51 = cos(t59) / 0.2e1;
t60 = pkin(7) - qJ(3);
t57 = cos(t60);
t47 = t57 / 0.2e1 + t51;
t50 = sin(t60) / 0.2e1;
t53 = sin(t59);
t44 = t53 / 0.2e1 + t50;
t121 = m(5) - t125;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t120 = -t122 * t68 + t126 * t71 + mrSges(5,1);
t104 = m(3) + m(4) + t121;
t119 = -t67 * mrSges(7,1) - t70 * mrSges(7,2) + t125 * pkin(9) + mrSges(5,2) - mrSges(6,3);
t69 = sin(qJ(3));
t116 = t30 * t69;
t115 = t31 * t69;
t114 = t43 * t69;
t40 = t44 * pkin(3);
t41 = t47 * pkin(3);
t99 = -pkin(3) * t115 + t40 * t107 - t74 * t41;
t98 = -pkin(3) * t114 + t66 * t40 + t80 * t41;
t83 = -pkin(3) * t116 - t40 * t113 - t73 * t41;
t79 = t97 / 0.2e1 + t90;
t78 = t96 / 0.2e1 + t89;
t76 = t62 * t78;
t72 = cos(qJ(3));
t65 = cos(pkin(7));
t61 = sin(pkin(7));
t52 = sin(t58);
t46 = t51 - t57 / 0.2e1;
t45 = t50 - t53 / 0.2e1;
t32 = -t80 * t61 + t66 * t65;
t26 = t65 * t107 + t74 * t61;
t25 = -t65 * t113 + t73 * t61;
t16 = t43 * t52 - t66 * t78 - t80 * t79;
t12 = -t112 * t76 + t31 * t52 + t74 * t79;
t9 = t30 * t52 + t64 * t76 + t73 * t79;
t6 = t17 * t71 + t32 * t68;
t4 = t13 * t71 + t26 * t68;
t2 = t10 * t71 + t25 * t68;
t1 = [(-m(2) - t104) * g(3) (-t66 * g(3) + (-g(1) * t112 + t64 * g(2)) * t62) * t104 (-(t66 * t44 + t80 * t47 - t114) * mrSges(4,1) - (-t43 * t72 + t80 * t45 + t66 * t46) * mrSges(4,2) - m(5) * t98 + t125 * (-t16 * pkin(4) + t98) + t119 * t17 + t120 * t16) * g(3) + (-(-t44 * t113 - t73 * t47 - t116) * mrSges(4,1) - (-t46 * t113 - t30 * t72 - t73 * t45) * mrSges(4,2) - m(5) * t83 + t125 * (-t9 * pkin(4) + t83) + t120 * t9 + t119 * t10) * g(2) + (-(t44 * t107 - t74 * t47 - t115) * mrSges(4,1) - (t46 * t107 - t31 * t72 - t74 * t45) * mrSges(4,2) - m(5) * t99 + t125 * (-t12 * pkin(4) + t99) + t119 * t13 + t120 * t12) * g(1), t121 * (-g(1) * t26 - g(2) * t25 - g(3) * t32) (t122 * t6 - t126 * (-t17 * t68 + t32 * t71)) * g(3) + (t122 * t2 - t126 * (-t10 * t68 + t25 * t71)) * g(2) + (t122 * t4 - t126 * (-t13 * t68 + t26 * t71)) * g(1), -g(1) * ((t12 * t70 - t4 * t67) * mrSges(7,1) + (-t12 * t67 - t4 * t70) * mrSges(7,2)) - g(2) * ((-t2 * t67 + t70 * t9) * mrSges(7,1) + (-t2 * t70 - t67 * t9) * mrSges(7,2)) - g(3) * ((t16 * t70 - t6 * t67) * mrSges(7,1) + (-t16 * t67 - t6 * t70) * mrSges(7,2))];
taug  = t1(:);
