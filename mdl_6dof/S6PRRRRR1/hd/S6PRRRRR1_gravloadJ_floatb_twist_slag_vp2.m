% Calculate Gravitation load on the joints for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:32:01
% EndTime: 2018-11-23 15:32:02
% DurationCPUTime: 0.84s
% Computational Cost: add. (1390->135), mult. (1264->173), div. (0->0), fcn. (1220->18), ass. (0->76)
t149 = mrSges(6,2) - mrSges(7,3);
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t148 = -t75 * mrSges(7,1) + t72 * mrSges(7,2) - mrSges(6,1);
t69 = sin(pkin(12));
t70 = sin(pkin(6));
t127 = t69 * t70;
t114 = cos(pkin(12));
t77 = cos(qJ(2));
t102 = t114 * t77;
t113 = pkin(6) - qJ(2);
t100 = sin(t113);
t112 = pkin(6) + qJ(2);
t62 = sin(t112);
t133 = t62 / 0.2e1;
t111 = t133 - t100 / 0.2e1;
t40 = -t111 * t69 + t102;
t68 = qJ(3) + qJ(4);
t63 = sin(t68);
t64 = cos(t68);
t147 = t64 * t127 - t40 * t63;
t101 = cos(t113);
t97 = cos(t112) / 0.2e1;
t51 = t97 - t101 / 0.2e1;
t71 = cos(pkin(6));
t145 = t51 * t63 + t64 * t71;
t144 = -m(6) - m(7);
t143 = -m(5) * pkin(3) - mrSges(4,1);
t65 = qJ(5) + t68;
t59 = sin(t65);
t60 = cos(t65);
t28 = t51 * t59 + t60 * t71;
t29 = -t51 * t60 + t59 * t71;
t142 = t148 * t28 + t149 * t29;
t13 = t127 * t60 - t40 * t59;
t14 = t127 * t59 + t40 * t60;
t141 = t148 * t13 + t149 * t14;
t103 = t70 * t114;
t126 = t69 * t77;
t37 = t111 * t114 + t126;
t11 = -t103 * t60 - t37 * t59;
t12 = -t103 * t59 + t37 * t60;
t140 = t148 * t11 + t149 * t12;
t139 = -t145 * mrSges(5,1) - (t51 * t64 - t63 * t71) * mrSges(5,2) + t142;
t138 = -t147 * mrSges(5,1) - (-t127 * t63 - t40 * t64) * mrSges(5,2) + t141;
t82 = -t103 * t64 - t37 * t63;
t137 = -t82 * mrSges(5,1) - (t103 * t63 - t37 * t64) * mrSges(5,2) + t140;
t76 = cos(qJ(3));
t66 = t76 * pkin(3);
t73 = sin(qJ(3));
t136 = mrSges(3,1) + m(5) * (t66 + pkin(2)) + t64 * mrSges(5,1) - t63 * mrSges(5,2) + m(4) * pkin(2) + t76 * mrSges(4,1) - t73 * mrSges(4,2) + (m(7) * pkin(5) - t148) * t60 + (m(7) * pkin(11) - t149) * t59;
t78 = -pkin(9) - pkin(8);
t135 = -m(4) * pkin(8) + m(5) * t78 - t72 * mrSges(7,1) - t75 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t52 = -pkin(3) * t73 - pkin(4) * t63;
t53 = pkin(4) * t64 + t66;
t119 = t53 * t127 + t40 * t52;
t115 = -t51 * t52 + t71 * t53;
t107 = t11 * pkin(5) + t12 * pkin(11);
t105 = t13 * pkin(5) + pkin(11) * t14;
t104 = t28 * pkin(5) + pkin(11) * t29;
t99 = t147 * pkin(4);
t98 = t145 * pkin(4);
t96 = t100 / 0.2e1;
t89 = t101 / 0.2e1 + t97;
t87 = -t103 * t53 + t37 * t52;
t86 = t96 - t62 / 0.2e1;
t80 = t82 * pkin(4);
t74 = sin(qJ(2));
t67 = -pkin(10) + t78;
t50 = t133 + t96;
t49 = pkin(2) + t53;
t41 = t69 * t86 + t102;
t39 = t114 * t74 + t69 * t89;
t38 = -t114 * t86 + t126;
t36 = -t114 * t89 + t69 * t74;
t1 = [(-m(2) - m(3) - m(4) - m(5) + t144) * g(3) (t144 * (t50 * t49 + t51 * t67) - t135 * t51 - t136 * t50) * g(3) + (t144 * (-t36 * t49 - t38 * t67) + t135 * t38 + t136 * t36) * g(2) + (t144 * (-t39 * t49 - t41 * t67) + t135 * t41 + t136 * t39) * g(1) (-(t51 * t76 - t71 * t73) * mrSges(4,2) - m(6) * t115 - m(7) * (t104 + t115) + t143 * (t51 * t73 + t71 * t76) + t139) * g(3) + (-(t103 * t73 - t37 * t76) * mrSges(4,2) - m(6) * t87 - m(7) * (t107 + t87) + t143 * (-t103 * t76 - t37 * t73) + t137) * g(2) + (-(-t127 * t73 - t40 * t76) * mrSges(4,2) - m(6) * t119 - m(7) * (t105 + t119) + t143 * (t127 * t76 - t40 * t73) + t138) * g(1) (-m(6) * t98 - m(7) * (t104 + t98) + t139) * g(3) + (-m(6) * t80 - m(7) * (t107 + t80) + t137) * g(2) + (-m(6) * t99 - m(7) * (t105 + t99) + t138) * g(1) (-m(7) * t104 + t142) * g(3) + (-m(7) * t107 + t140) * g(2) + (-m(7) * t105 + t141) * g(1), -g(1) * ((-t14 * t72 + t39 * t75) * mrSges(7,1) + (-t14 * t75 - t39 * t72) * mrSges(7,2)) - g(2) * ((-t12 * t72 + t36 * t75) * mrSges(7,1) + (-t12 * t75 - t36 * t72) * mrSges(7,2)) - g(3) * ((-t29 * t72 - t50 * t75) * mrSges(7,1) + (-t29 * t75 + t50 * t72) * mrSges(7,2))];
taug  = t1(:);
