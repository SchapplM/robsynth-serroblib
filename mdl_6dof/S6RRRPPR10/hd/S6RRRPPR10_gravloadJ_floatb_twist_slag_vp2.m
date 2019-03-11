% Calculate Gravitation load on the joints for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:29
% EndTime: 2019-03-09 16:21:33
% DurationCPUTime: 1.47s
% Computational Cost: add. (660->121), mult. (1502->155), div. (0->0), fcn. (1760->12), ass. (0->67)
t44 = sin(pkin(11));
t46 = cos(pkin(11));
t127 = (-m(7) * pkin(5) - mrSges(6,1)) * t44 - t46 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t43 = pkin(11) + qJ(6);
t40 = sin(t43);
t41 = cos(t43);
t63 = t40 * mrSges(7,1) + t41 * mrSges(7,2);
t134 = -t63 + t127;
t133 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) + m(7) * (-pkin(10) - qJ(5));
t132 = -mrSges(6,3) - mrSges(7,3);
t112 = m(6) + m(7);
t108 = m(5) + t112;
t131 = qJ(4) * t108;
t57 = t132 + t133;
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t123 = -pkin(3) * t51 - qJ(4) * t48;
t119 = pkin(3) * t108 - t57;
t64 = t41 * mrSges(7,1) - t40 * mrSges(7,2);
t116 = -t46 * mrSges(6,1) + t44 * mrSges(6,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t39 = pkin(5) * t46 + pkin(4);
t115 = -m(6) * (pkin(4) + pkin(9)) + t116 - m(7) * (pkin(9) + t39) - t64;
t114 = -t134 * t48 - t51 * t57 + mrSges(3,1);
t113 = -m(6) * pkin(4) - m(7) * t39 + t116;
t49 = sin(qJ(2));
t50 = sin(qJ(1));
t86 = cos(pkin(6));
t98 = cos(qJ(2));
t66 = t86 * t98;
t99 = cos(qJ(1));
t24 = t49 * t50 - t99 * t66;
t110 = t123 * t24;
t26 = t99 * t49 + t50 * t66;
t109 = t123 * t26;
t107 = -(m(4) + t108) * pkin(9) + t113;
t104 = t134 - t131;
t45 = sin(pkin(6));
t96 = t45 * t49;
t95 = t45 * t50;
t94 = t45 * t51;
t79 = t45 * t98;
t90 = pkin(2) * t79 + pkin(9) * t96;
t89 = t99 * pkin(1) + pkin(8) * t95;
t72 = t49 * t86;
t27 = -t50 * t72 + t99 * t98;
t82 = t27 * pkin(2) + t89;
t80 = t45 * t99;
t78 = t48 * t98;
t77 = t51 * t98;
t75 = -pkin(1) * t50 + pkin(8) * t80;
t18 = t24 * pkin(2);
t25 = t50 * t98 + t99 * t72;
t74 = t25 * pkin(9) - t18;
t20 = t26 * pkin(2);
t73 = t27 * pkin(9) - t20;
t8 = t25 * t51 - t48 * t80;
t69 = t45 * t77;
t67 = -t25 * pkin(2) + t75;
t7 = t25 * t48 + t51 * t80;
t53 = t127 - t131;
t23 = t86 * t48 + t49 * t94;
t22 = t48 * t96 - t86 * t51;
t12 = t27 * t51 + t48 * t95;
t11 = t27 * t48 - t50 * t94;
t2 = t11 * t40 + t26 * t41;
t1 = t11 * t41 - t26 * t40;
t3 = [(-t99 * mrSges(2,1) - m(3) * t89 - t27 * mrSges(3,1) - m(4) * t82 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t45 + mrSges(2,2)) * t50 + t57 * t12 + t53 * t11 + t107 * t26 - t108 * (t12 * pkin(3) + t82)) * g(2) + (t50 * mrSges(2,1) + t99 * mrSges(2,2) - m(3) * t75 + t25 * mrSges(3,1) - mrSges(3,3) * t80 - m(4) * t67 - t57 * t8 - (t53 - t63) * t7 + (-t107 + t64) * t24 + t108 * (pkin(3) * t8 - t67)) * g(1) (-m(4) * t90 + t132 * t69 - t108 * (t45 * qJ(4) * t78 + pkin(3) * t69 + t90) + (-mrSges(3,1) * t98 + t133 * t77 + (t113 - t64) * t49 + t134 * t78) * t45) * g(3) + (-m(4) * t74 - m(5) * (t74 + t110) - t112 * (-t18 + t110) + t115 * t25 + t114 * t24) * g(2) + (-m(4) * t73 - m(5) * (t73 + t109) - t112 * (-t20 + t109) + t115 * t27 + t114 * t26) * g(1) (t104 * t23 + t119 * t22) * g(3) + (t104 * t8 + t119 * t7) * g(2) + (t104 * t12 + t119 * t11) * g(1), t108 * (-g(1) * t11 - g(2) * t7 - g(3) * t22) t112 * (-g(1) * t12 - g(2) * t8 - g(3) * t23) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t24 * t40 + t41 * t7) * mrSges(7,1) + (-t24 * t41 - t40 * t7) * mrSges(7,2)) - g(3) * ((t22 * t41 + t40 * t79) * mrSges(7,1) + (-t22 * t40 + t41 * t79) * mrSges(7,2))];
taug  = t3(:);
