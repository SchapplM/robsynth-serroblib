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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:13
% EndTime: 2019-03-08 23:00:16
% DurationCPUTime: 0.93s
% Computational Cost: add. (702->123), mult. (963->168), div. (0->0), fcn. (1077->14), ass. (0->61)
t135 = mrSges(6,2) - mrSges(7,3);
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t134 = -mrSges(7,1) * t65 + mrSges(7,2) * t62 - mrSges(6,1);
t100 = cos(pkin(6));
t61 = sin(pkin(6));
t64 = sin(qJ(2));
t113 = t61 * t64;
t59 = qJ(3) + qJ(4);
t55 = sin(t59);
t56 = cos(t59);
t133 = t100 * t56 - t113 * t55;
t60 = sin(pkin(11));
t114 = t60 * t61;
t67 = cos(qJ(2));
t90 = t60 * t100;
t99 = cos(pkin(11));
t40 = -t64 * t90 + t67 * t99;
t132 = t114 * t56 - t40 * t55;
t68 = -pkin(9) - pkin(8);
t120 = -m(4) * pkin(8) + m(5) * t68 - t62 * mrSges(7,1) - t65 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t54 = pkin(12) + t59;
t50 = sin(t54);
t51 = cos(t54);
t66 = cos(qJ(3));
t57 = t66 * pkin(3);
t63 = sin(qJ(3));
t130 = -mrSges(3,1) - m(5) * (t57 + pkin(2)) - t56 * mrSges(5,1) + t55 * mrSges(5,2) - m(4) * pkin(2) - t66 * mrSges(4,1) + t63 * mrSges(4,2) + (-m(7) * pkin(5) + t134) * t51 + (-m(7) * pkin(10) + t135) * t50;
t129 = m(6) + m(7);
t127 = -m(5) * pkin(3) - mrSges(4,1);
t28 = t100 * t51 - t113 * t50;
t29 = t100 * t50 + t113 * t51;
t124 = -t133 * mrSges(5,1) - (-t100 * t55 - t113 * t56) * mrSges(5,2) + t135 * t29 + t134 * t28;
t13 = t114 * t51 - t40 * t50;
t14 = t114 * t50 + t40 * t51;
t123 = -t132 * mrSges(5,1) - (-t114 * t55 - t40 * t56) * mrSges(5,2) + t135 * t14 + t134 * t13;
t81 = t100 * t99;
t38 = t60 * t67 + t64 * t81;
t89 = t61 * t99;
t11 = -t38 * t50 - t51 * t89;
t12 = t38 * t51 - t50 * t89;
t74 = -t38 * t55 - t56 * t89;
t122 = -t74 * mrSges(5,1) - (-t38 * t56 + t55 * t89) * mrSges(5,2) + t135 * t12 + t134 * t11;
t112 = t61 * t66;
t111 = t61 * t67;
t45 = -pkin(3) * t63 - pkin(4) * t55;
t46 = pkin(4) * t56 + t57;
t104 = t114 * t46 + t40 * t45;
t101 = t100 * t46 + t113 * t45;
t94 = pkin(5) * t11 + t12 * pkin(10);
t92 = pkin(5) * t13 + pkin(10) * t14;
t91 = pkin(5) * t28 + pkin(10) * t29;
t87 = t132 * pkin(4);
t82 = t133 * pkin(4);
t78 = t38 * t45 - t46 * t89;
t70 = t74 * pkin(4);
t58 = -qJ(5) + t68;
t44 = pkin(2) + t46;
t39 = t64 * t99 + t67 * t90;
t37 = t60 * t64 - t67 * t81;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t129) * g(3) (-t129 * (-t37 * t44 - t38 * t58) + t120 * t38 - t130 * t37) * g(2) + (-t129 * (-t39 * t44 - t40 * t58) + t120 * t40 - t130 * t39) * g(1) + (-t129 * t44 * t111 + (t130 * t67 + (t129 * t58 + t120) * t64) * t61) * g(3) (-(-t100 * t63 - t112 * t64) * mrSges(4,2) - m(6) * t101 - m(7) * (t91 + t101) + t127 * (t100 * t66 - t113 * t63) + t124) * g(3) + (-(-t38 * t66 + t63 * t89) * mrSges(4,2) - m(6) * t78 - m(7) * (t78 + t94) + t127 * (-t38 * t63 - t66 * t89) + t122) * g(2) + (-(-t114 * t63 - t40 * t66) * mrSges(4,2) - m(6) * t104 - m(7) * (t92 + t104) + t127 * (t112 * t60 - t40 * t63) + t123) * g(1) (-m(6) * t82 - m(7) * (t82 + t91) + t124) * g(3) + (-m(6) * t70 - m(7) * (t70 + t94) + t122) * g(2) + (-m(6) * t87 - m(7) * (t87 + t92) + t123) * g(1), t129 * (-g(1) * t39 - g(2) * t37 + g(3) * t111) -g(1) * ((-t14 * t62 + t39 * t65) * mrSges(7,1) + (-t14 * t65 - t39 * t62) * mrSges(7,2)) - g(2) * ((-t12 * t62 + t37 * t65) * mrSges(7,1) + (-t12 * t65 - t37 * t62) * mrSges(7,2)) - g(3) * ((-t111 * t65 - t29 * t62) * mrSges(7,1) + (t111 * t62 - t29 * t65) * mrSges(7,2))];
taug  = t1(:);
