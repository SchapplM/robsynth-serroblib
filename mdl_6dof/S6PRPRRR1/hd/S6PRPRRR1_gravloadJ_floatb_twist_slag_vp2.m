% Calculate Gravitation load on the joints for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:03:25
% EndTime: 2018-11-23 15:03:26
% DurationCPUTime: 0.87s
% Computational Cost: add. (1231->122), mult. (972->165), div. (0->0), fcn. (897->22), ass. (0->77)
t138 = mrSges(6,2) - mrSges(7,3);
t70 = sin(qJ(6));
t73 = cos(qJ(6));
t137 = -t73 * mrSges(7,1) + t70 * mrSges(7,2) - mrSges(6,1);
t64 = pkin(6) - qJ(2);
t51 = sin(t64) / 0.2e1;
t63 = pkin(6) + qJ(2);
t55 = sin(t63);
t136 = t55 / 0.2e1 + t51;
t67 = sin(pkin(6));
t74 = cos(qJ(4));
t115 = t67 * t74;
t66 = sin(pkin(11));
t71 = sin(qJ(4));
t62 = qJ(2) + pkin(12);
t100 = pkin(6) + t62;
t83 = sin(t100) / 0.2e1;
t101 = pkin(6) - t62;
t89 = sin(t101);
t37 = t83 - t89 / 0.2e1;
t57 = cos(t62);
t68 = cos(pkin(11));
t93 = -t37 * t66 + t57 * t68;
t135 = t66 * t115 - t93 * t71;
t84 = cos(t101) / 0.2e1;
t90 = cos(t100);
t39 = t84 - t90 / 0.2e1;
t69 = cos(pkin(6));
t133 = -t39 * t71 + t69 * t74;
t132 = -m(4) - m(5);
t131 = -m(6) - m(7);
t65 = qJ(4) + qJ(5);
t60 = sin(t65);
t61 = cos(t65);
t24 = -t39 * t60 + t61 * t69;
t25 = t39 * t61 + t60 * t69;
t130 = t137 * t24 + t138 * t25;
t119 = t66 * t67;
t13 = t119 * t61 - t60 * t93;
t14 = t119 * t60 + t61 * t93;
t129 = t137 * t13 + t138 * t14;
t117 = t67 * t68;
t94 = t68 * t37 + t57 * t66;
t11 = -t117 * t61 - t60 * t94;
t12 = -t117 * t60 + t61 * t94;
t128 = t137 * t11 + t138 * t12;
t52 = cos(t63) / 0.2e1;
t59 = cos(t64);
t42 = t59 / 0.2e1 + t52;
t127 = m(5) * pkin(3) + t74 * mrSges(5,1) - t71 * mrSges(5,2) + mrSges(4,1) + (m(7) * pkin(5) - t137) * t61 + (m(7) * pkin(10) - t138) * t60;
t126 = -m(5) * pkin(8) - t70 * mrSges(7,1) - t73 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t72 = sin(qJ(2));
t118 = t66 * t72;
t116 = t67 * t71;
t114 = t68 * t72;
t109 = t136 * pkin(2);
t107 = -t131 - t132;
t105 = t11 * pkin(5) + pkin(10) * t12;
t104 = t13 * pkin(5) + pkin(10) * t14;
t103 = t24 * pkin(5) + pkin(10) * t25;
t40 = t42 * pkin(2);
t99 = -pkin(2) * t118 + t68 * t40;
t98 = t135 * pkin(4);
t97 = t133 * pkin(4);
t87 = -pkin(2) * t114 - t40 * t66;
t86 = -t115 * t68 - t71 * t94;
t80 = t86 * pkin(4);
t78 = t90 / 0.2e1 + t84;
t76 = -pkin(9) - pkin(8);
t75 = cos(qJ(2));
t54 = sin(t62);
t53 = pkin(4) * t74 + pkin(3);
t41 = t51 - t55 / 0.2e1;
t38 = t89 / 0.2e1 + t83;
t29 = t68 * t54 + t66 * t78;
t26 = t54 * t66 - t68 * t78;
t1 = [(-m(2) - m(3) - t107) * g(3) (-t136 * mrSges(3,1) - (t52 - t59 / 0.2e1) * mrSges(3,2) + t132 * t109 + t131 * (t38 * t53 - t39 * t76 + t109) + t126 * t39 - t127 * t38) * g(3) + (-(t42 * t68 - t118) * mrSges(3,1) - (t41 * t68 - t66 * t75) * mrSges(3,2) + t132 * t99 + t131 * (-t26 * t53 - t76 * t94 + t99) + t126 * t94 + t127 * t26) * g(2) + (-(-t42 * t66 - t114) * mrSges(3,1) - (-t41 * t66 - t68 * t75) * mrSges(3,2) + t132 * t87 + t131 * (-t29 * t53 - t76 * t93 + t87) + t126 * t93 + t127 * t29) * g(1) (-g(3) * t69 + (-g(1) * t66 + g(2) * t68) * t67) * t107 (-t133 * mrSges(5,1) - (-t39 * t74 - t69 * t71) * mrSges(5,2) - m(6) * t97 - m(7) * (t103 + t97) + t130) * g(3) + (-t86 * mrSges(5,1) - (t116 * t68 - t74 * t94) * mrSges(5,2) - m(6) * t80 - m(7) * (t105 + t80) + t128) * g(2) + (-t135 * mrSges(5,1) - (-t116 * t66 - t74 * t93) * mrSges(5,2) - m(6) * t98 - m(7) * (t104 + t98) + t129) * g(1) (-m(7) * t103 + t130) * g(3) + (-m(7) * t105 + t128) * g(2) + (-m(7) * t104 + t129) * g(1), -g(1) * ((-t14 * t70 + t29 * t73) * mrSges(7,1) + (-t14 * t73 - t29 * t70) * mrSges(7,2)) - g(2) * ((-t12 * t70 + t26 * t73) * mrSges(7,1) + (-t12 * t73 - t26 * t70) * mrSges(7,2)) - g(3) * ((-t25 * t70 - t38 * t73) * mrSges(7,1) + (-t25 * t73 + t38 * t70) * mrSges(7,2))];
taug  = t1(:);
