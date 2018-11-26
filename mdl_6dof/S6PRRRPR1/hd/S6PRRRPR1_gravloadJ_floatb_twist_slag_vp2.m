% Calculate Gravitation load on the joints for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:22:19
% EndTime: 2018-11-23 15:22:20
% DurationCPUTime: 0.83s
% Computational Cost: add. (1215->132), mult. (1134->171), div. (0->0), fcn. (1077->18), ass. (0->73)
t148 = mrSges(6,2) - mrSges(7,3);
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t147 = -t75 * mrSges(7,1) + t72 * mrSges(7,2) - mrSges(6,1);
t69 = sin(pkin(11));
t70 = sin(pkin(6));
t130 = t69 * t70;
t117 = cos(pkin(11));
t77 = cos(qJ(2));
t106 = t117 * t77;
t116 = pkin(6) - qJ(2);
t104 = sin(t116);
t115 = pkin(6) + qJ(2);
t103 = sin(t115);
t98 = t103 / 0.2e1;
t89 = t98 - t104 / 0.2e1;
t40 = -t69 * t89 + t106;
t68 = qJ(3) + qJ(4);
t64 = sin(t68);
t65 = cos(t68);
t146 = t65 * t130 - t40 * t64;
t100 = cos(t115) / 0.2e1;
t105 = cos(t116);
t51 = t100 - t105 / 0.2e1;
t71 = cos(pkin(6));
t144 = t51 * t64 + t65 * t71;
t143 = m(6) + m(7);
t142 = -m(5) * pkin(3) - mrSges(4,1);
t63 = pkin(12) + t68;
t59 = sin(t63);
t60 = cos(t63);
t28 = t51 * t59 + t60 * t71;
t29 = -t51 * t60 + t59 * t71;
t141 = -t144 * mrSges(5,1) - (t51 * t65 - t64 * t71) * mrSges(5,2) + t148 * t29 + t147 * t28;
t13 = t130 * t60 - t40 * t59;
t14 = t130 * t59 + t40 * t60;
t140 = -t146 * mrSges(5,1) - (-t130 * t64 - t40 * t65) * mrSges(5,2) + t148 * t14 + t147 * t13;
t107 = t70 * t117;
t129 = t69 * t77;
t37 = t117 * t89 + t129;
t11 = -t107 * t60 - t37 * t59;
t12 = -t107 * t59 + t37 * t60;
t82 = -t107 * t65 - t37 * t64;
t139 = -t82 * mrSges(5,1) - (t107 * t64 - t37 * t65) * mrSges(5,2) + t148 * t12 + t147 * t11;
t76 = cos(qJ(3));
t66 = t76 * pkin(3);
t73 = sin(qJ(3));
t138 = mrSges(3,1) + m(5) * (t66 + pkin(2)) + t65 * mrSges(5,1) - t64 * mrSges(5,2) + m(4) * pkin(2) + t76 * mrSges(4,1) - t73 * mrSges(4,2) + (m(7) * pkin(5) - t147) * t60 + (m(7) * pkin(10) - t148) * t59;
t78 = -pkin(9) - pkin(8);
t137 = -m(4) * pkin(8) + m(5) * t78 - t72 * mrSges(7,1) - t75 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t52 = -pkin(3) * t73 - pkin(4) * t64;
t53 = pkin(4) * t65 + t66;
t122 = t53 * t130 + t40 * t52;
t118 = -t51 * t52 + t71 * t53;
t111 = t11 * pkin(5) + t12 * pkin(10);
t109 = t13 * pkin(5) + pkin(10) * t14;
t108 = t28 * pkin(5) + pkin(10) * t29;
t102 = t146 * pkin(4);
t101 = t144 * pkin(4);
t99 = t104 / 0.2e1;
t91 = t105 / 0.2e1 + t100;
t90 = t99 - t103 / 0.2e1;
t87 = -t107 * t53 + t37 * t52;
t80 = t82 * pkin(4);
t74 = sin(qJ(2));
t67 = -qJ(5) + t78;
t50 = t98 + t99;
t49 = pkin(2) + t53;
t41 = t69 * t90 + t106;
t39 = t117 * t74 + t69 * t91;
t38 = -t117 * t90 + t129;
t36 = -t117 * t91 + t69 * t74;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t143) * g(3) (-t143 * (t50 * t49 + t51 * t67) - t137 * t51 - t138 * t50) * g(3) + (-t143 * (-t36 * t49 - t38 * t67) + t137 * t38 + t138 * t36) * g(2) + (-t143 * (-t39 * t49 - t41 * t67) + t137 * t41 + t138 * t39) * g(1) (-(t51 * t76 - t71 * t73) * mrSges(4,2) - m(6) * t118 - m(7) * (t108 + t118) + t142 * (t51 * t73 + t71 * t76) + t141) * g(3) + (-(t107 * t73 - t37 * t76) * mrSges(4,2) - m(6) * t87 - m(7) * (t111 + t87) + t142 * (-t107 * t76 - t37 * t73) + t139) * g(2) + (-(-t130 * t73 - t40 * t76) * mrSges(4,2) - m(6) * t122 - m(7) * (t109 + t122) + t142 * (t130 * t76 - t40 * t73) + t140) * g(1) (-m(6) * t101 - m(7) * (t101 + t108) + t141) * g(3) + (-m(6) * t80 - m(7) * (t111 + t80) + t139) * g(2) + (-m(6) * t102 - m(7) * (t102 + t109) + t140) * g(1), t143 * (-g(1) * t39 - g(2) * t36 + g(3) * t50) -g(1) * ((-t14 * t72 + t39 * t75) * mrSges(7,1) + (-t14 * t75 - t39 * t72) * mrSges(7,2)) - g(2) * ((-t12 * t72 + t36 * t75) * mrSges(7,1) + (-t12 * t75 - t36 * t72) * mrSges(7,2)) - g(3) * ((-t29 * t72 - t50 * t75) * mrSges(7,1) + (-t29 * t75 + t50 * t72) * mrSges(7,2))];
taug  = t1(:);
