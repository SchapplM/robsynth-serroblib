% Calculate Gravitation load on the joints for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:22
% EndTime: 2019-03-09 14:08:25
% DurationCPUTime: 1.23s
% Computational Cost: add. (868->144), mult. (1194->190), div. (0->0), fcn. (1354->14), ass. (0->68)
t139 = mrSges(6,2) - mrSges(7,3);
t69 = sin(qJ(6));
t72 = cos(qJ(6));
t138 = t72 * mrSges(7,1) - t69 * mrSges(7,2) + mrSges(6,1);
t140 = -m(7) * pkin(5) - t138;
t92 = -m(7) * pkin(11) + t139;
t68 = -pkin(9) - qJ(3);
t76 = -m(4) * qJ(3) + m(5) * t68 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t130 = t69 * mrSges(7,1) + t72 * mrSges(7,2) - t76;
t107 = cos(pkin(6));
t66 = sin(pkin(6));
t70 = sin(qJ(2));
t116 = t66 * t70;
t64 = pkin(12) + qJ(4);
t58 = sin(t64);
t59 = cos(t64);
t137 = t107 * t59 - t58 * t116;
t71 = sin(qJ(1));
t115 = t66 * t71;
t118 = cos(qJ(1));
t73 = cos(qJ(2));
t94 = t71 * t107;
t42 = t118 * t73 - t70 * t94;
t19 = t59 * t115 - t42 * t58;
t67 = cos(pkin(12));
t57 = t67 * pkin(3) + pkin(2);
t65 = sin(pkin(12));
t135 = -m(4) * pkin(2) - m(5) * t57 - t67 * mrSges(4,1) + t65 * mrSges(4,2) - mrSges(3,1);
t60 = qJ(5) + t64;
t55 = sin(t60);
t56 = cos(t60);
t134 = -t59 * mrSges(5,1) + t58 * mrSges(5,2) + t140 * t56 + t92 * t55 + t135;
t133 = -m(4) - m(5);
t132 = m(6) + m(7);
t30 = t107 * t56 - t116 * t55;
t31 = t107 * t55 + t116 * t56;
t129 = -t138 * t30 + t139 * t31;
t17 = -t115 * t56 + t42 * t55;
t18 = t115 * t55 + t42 * t56;
t128 = t138 * t17 + t139 * t18;
t102 = t66 * t118;
t88 = t107 * t118;
t40 = t70 * t88 + t71 * t73;
t13 = -t56 * t102 - t40 * t55;
t14 = -t55 * t102 + t40 * t56;
t127 = -t138 * t13 + t139 * t14;
t114 = t66 * t73;
t108 = t118 * pkin(1) + pkin(8) * t115;
t100 = -t71 * pkin(1) + pkin(8) * t102;
t99 = t13 * pkin(5) + t14 * pkin(11);
t98 = -t17 * pkin(5) + pkin(11) * t18;
t97 = t30 * pkin(5) + pkin(11) * t31;
t95 = t58 * t102 - t40 * t59;
t91 = t65 * t102;
t90 = t19 * pkin(4);
t41 = t118 * t70 + t73 * t94;
t43 = pkin(4) * t59 + t57;
t48 = pkin(3) * t65 + pkin(4) * t58;
t63 = -pkin(10) + t68;
t89 = t48 * t115 - t41 * t63 + t42 * t43 + t108;
t83 = t137 * pkin(4);
t78 = t102 * t59 + t40 * t58;
t75 = t78 * pkin(4);
t39 = t70 * t71 - t73 * t88;
t20 = t115 * t58 + t42 * t59;
t2 = t18 * t72 + t41 * t69;
t1 = -t18 * t69 + t41 * t72;
t3 = [(t71 * mrSges(2,1) + t118 * mrSges(2,2) - m(3) * t100 + t40 * mrSges(3,1) - mrSges(3,3) * t102 - m(4) * (-pkin(2) * t40 + t100) - (-t40 * t67 + t91) * mrSges(4,1) - (t102 * t67 + t40 * t65) * mrSges(4,2) - m(5) * (pkin(3) * t91 - t40 * t57 + t100) - t95 * mrSges(5,1) - t78 * mrSges(5,2) + t92 * t13 - t140 * t14 + t130 * t39 + t132 * (-t48 * t102 - t39 * t63 + t40 * t43 - t100)) * g(1) + (-t118 * mrSges(2,1) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t89 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t89) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t135 * t42 + (mrSges(2,2) + (-mrSges(4,2) * t67 - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t65) * t66) * t71 + t92 * t17 + t76 * t41 + (-m(3) + t133) * t108) * g(2) (-t132 * (-t39 * t43 - t40 * t63) - t130 * t40 - t134 * t39) * g(2) + (-t132 * (-t41 * t43 - t42 * t63) - t130 * t42 - t134 * t41) * g(1) + (-t132 * t43 * t114 + (t134 * t73 + (t132 * t63 - t130) * t70) * t66) * g(3) (-g(1) * t41 - g(2) * t39 + g(3) * t114) * (t132 - t133) (-t137 * mrSges(5,1) - (-t107 * t58 - t116 * t59) * mrSges(5,2) - m(6) * t83 - m(7) * (t83 + t97) + t129) * g(3) + (t78 * mrSges(5,1) - t95 * mrSges(5,2) + m(6) * t75 - m(7) * (-t75 + t99) + t127) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t90 - m(7) * (t90 + t98) + t128) * g(1) (-m(7) * t97 + t129) * g(3) + (-m(7) * t99 + t127) * g(2) + (-m(7) * t98 + t128) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t69 + t39 * t72) * mrSges(7,1) + (-t14 * t72 - t39 * t69) * mrSges(7,2)) - g(3) * ((-t114 * t72 - t31 * t69) * mrSges(7,1) + (t114 * t69 - t31 * t72) * mrSges(7,2))];
taug  = t3(:);
