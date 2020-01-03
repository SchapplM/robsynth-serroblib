% Calculate Gravitation load on the joints for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:19
% EndTime: 2019-12-31 21:43:22
% DurationCPUTime: 0.89s
% Computational Cost: add. (404->107), mult. (983->150), div. (0->0), fcn. (1144->10), ass. (0->63)
t41 = sin(qJ(5));
t45 = cos(qJ(5));
t106 = t41 * mrSges(6,1) + t45 * mrSges(6,2);
t100 = mrSges(4,2) - mrSges(5,3);
t101 = m(5) + m(6);
t50 = -qJ(4) * t101 + t100;
t90 = t50 - t106;
t107 = -mrSges(4,1) + mrSges(5,2);
t94 = -m(6) * pkin(9) - mrSges(6,3);
t91 = -t94 - t107;
t103 = pkin(3) * t101 + t91;
t102 = m(6) * pkin(4);
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t47 = cos(qJ(2));
t68 = cos(pkin(5));
t84 = cos(qJ(1));
t57 = t68 * t84;
t24 = t43 * t44 - t47 * t57;
t42 = sin(qJ(3));
t69 = qJ(4) * t42;
t46 = cos(qJ(3));
t83 = t24 * t46;
t99 = -pkin(3) * t83 - t24 * t69;
t61 = t44 * t68;
t26 = t43 * t84 + t47 * t61;
t82 = t26 * t46;
t98 = -pkin(3) * t82 - t26 * t69;
t97 = t100 * t42 + t107 * t46 - mrSges(3,1);
t96 = mrSges(3,2) - mrSges(5,1) - mrSges(4,3);
t93 = -pkin(8) * (m(4) + t101) - t102 + t96;
t92 = t45 * mrSges(6,1) - t41 * mrSges(6,2);
t89 = -t92 + t96;
t88 = t46 * mrSges(6,3) + t106 * t42 - t97;
t87 = -m(6) * (pkin(4) + pkin(8)) + t89;
t40 = sin(pkin(5));
t81 = t40 * t43;
t80 = t40 * t44;
t79 = t40 * t46;
t78 = t40 * t47;
t76 = t41 * t47;
t74 = t45 * t47;
t72 = t46 * t47;
t71 = pkin(2) * t78 + pkin(8) * t81;
t70 = t84 * pkin(1) + pkin(7) * t80;
t27 = -t43 * t61 + t47 * t84;
t66 = t27 * pkin(2) + t70;
t65 = t40 * t84;
t64 = -pkin(1) * t44 + pkin(7) * t65;
t18 = t24 * pkin(2);
t25 = t43 * t57 + t44 * t47;
t63 = pkin(8) * t25 - t18;
t20 = t26 * pkin(2);
t62 = pkin(8) * t27 - t20;
t8 = t25 * t46 - t42 * t65;
t58 = -t25 * pkin(2) + t64;
t7 = t25 * t42 + t46 * t65;
t22 = t42 * t81 - t46 * t68;
t12 = t27 * t46 + t42 * t80;
t11 = t27 * t42 - t44 * t79;
t2 = t11 * t41 + t26 * t45;
t1 = t11 * t45 - t26 * t41;
t3 = [(-t84 * mrSges(2,1) - m(3) * t70 - t27 * mrSges(3,1) - m(4) * t66 - t2 * mrSges(6,1) - t1 * mrSges(6,2) + (-mrSges(3,3) * t40 + mrSges(2,2)) * t44 + t50 * t11 - t91 * t12 + t93 * t26 - t101 * (t12 * pkin(3) + t66)) * g(2) + (t44 * mrSges(2,1) + t84 * mrSges(2,2) - m(3) * t64 + t25 * mrSges(3,1) - mrSges(3,3) * t65 - m(4) * t58 + t91 * t8 - t90 * t7 + (t92 - t93) * t24 + t101 * (pkin(3) * t8 - t58)) * g(1), (-m(4) * t71 - t101 * (t40 * pkin(3) * t72 + t69 * t78 + t71) + (t94 * t72 + (-t76 * mrSges(6,1) - t74 * mrSges(6,2)) * t42 + t97 * t47 + (t89 - t102) * t43) * t40) * g(3) + (-m(4) * t63 - m(5) * (t63 + t99) - m(6) * (-pkin(9) * t83 - t18 + t99) + t87 * t25 + t88 * t24) * g(2) + (-m(4) * t62 - m(5) * (t62 + t98) - m(6) * (-pkin(9) * t82 - t20 + t98) + t87 * t27 + t88 * t26) * g(1), (t90 * (t42 * t68 + t43 * t79) + t103 * t22) * g(3) + (t103 * t7 + t90 * t8) * g(2) + (t103 * t11 + t90 * t12) * g(1), t101 * (-g(1) * t11 - g(2) * t7 - g(3) * t22), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((-t24 * t41 + t45 * t7) * mrSges(6,1) + (-t24 * t45 - t41 * t7) * mrSges(6,2)) - g(3) * ((t22 * t45 + t40 * t76) * mrSges(6,1) + (-t22 * t41 + t40 * t74) * mrSges(6,2))];
taug = t3(:);
