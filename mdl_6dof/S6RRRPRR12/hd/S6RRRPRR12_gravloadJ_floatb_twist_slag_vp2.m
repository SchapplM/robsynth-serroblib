% Calculate Gravitation load on the joints for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:12
% EndTime: 2019-03-09 19:38:17
% DurationCPUTime: 1.67s
% Computational Cost: add. (839->128), mult. (1597->164), div. (0->0), fcn. (1881->14), ass. (0->65)
t51 = -pkin(10) - qJ(4);
t125 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(7,3) + m(7) * (-pkin(11) + t51) + m(6) * t51 - mrSges(6,3);
t50 = cos(pkin(12));
t40 = t50 * pkin(4) + pkin(3);
t47 = pkin(12) + qJ(5);
t42 = cos(t47);
t48 = sin(pkin(12));
t103 = -m(5) * pkin(3) - t50 * mrSges(5,1) + t48 * mrSges(5,2) - mrSges(4,1) - m(7) * (pkin(5) * t42 + t40) - m(6) * t40;
t43 = qJ(6) + t47;
t38 = sin(t43);
t39 = cos(t43);
t41 = sin(t47);
t136 = -t42 * mrSges(6,1) - t39 * mrSges(7,1) + t41 * mrSges(6,2) + t38 * mrSges(7,2) + t103;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t143 = t125 * t52 + t136 * t54 - mrSges(3,1);
t119 = -m(4) - m(5);
t138 = -t41 * mrSges(6,1) - t38 * mrSges(7,1) - t42 * mrSges(6,2) - t39 * mrSges(7,2);
t137 = t50 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3);
t118 = -m(6) - m(7);
t126 = -t118 - t119;
t135 = pkin(2) * t126 - t143;
t99 = pkin(4) * t48;
t28 = pkin(5) * t41 + t99;
t104 = t48 * mrSges(5,1) + m(7) * (pkin(9) + t28) + m(6) * (pkin(9) + t99) + t137;
t120 = t119 * pkin(9) - t104 + t138;
t109 = -m(7) * pkin(5) - mrSges(6,1);
t53 = sin(qJ(2));
t85 = cos(pkin(6));
t96 = cos(qJ(1));
t72 = t85 * t96;
t94 = sin(qJ(1));
t95 = cos(qJ(2));
t24 = t53 * t72 + t94 * t95;
t49 = sin(pkin(6));
t82 = t49 * t96;
t12 = t24 * t54 - t52 * t82;
t23 = t53 * t94 - t72 * t95;
t92 = t23 * t39;
t93 = t23 * t38;
t101 = (-t12 * t38 + t92) * mrSges(7,1) + (-t12 * t39 - t93) * mrSges(7,2);
t71 = t85 * t94;
t26 = -t53 * t71 + t95 * t96;
t80 = t49 * t94;
t16 = t26 * t54 + t52 * t80;
t25 = t53 * t96 + t71 * t95;
t5 = -t16 * t38 + t25 * t39;
t6 = t16 * t39 + t25 * t38;
t100 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t88 = t49 * t53;
t22 = t52 * t85 + t54 * t88;
t81 = t49 * t95;
t97 = (-t22 * t38 - t39 * t81) * mrSges(7,1) + (-t22 * t39 + t38 * t81) * mrSges(7,2);
t91 = t23 * t41;
t90 = t23 * t42;
t86 = t96 * pkin(1) + pkin(8) * t80;
t84 = t26 * pkin(2) + t86;
t73 = -pkin(1) * t94 + pkin(8) * t82;
t7 = -t16 * t41 + t25 * t42;
t67 = -t24 * pkin(2) + t73;
t11 = t24 * t52 + t54 * t82;
t21 = t52 * t88 - t54 * t85;
t15 = t26 * t52 - t54 * t80;
t8 = t16 * t42 + t25 * t41;
t1 = [(-m(3) * t86 - t96 * mrSges(2,1) - t26 * mrSges(3,1) - t8 * mrSges(6,1) - t6 * mrSges(7,1) + t94 * mrSges(2,2) - t7 * mrSges(6,2) - t5 * mrSges(7,2) - mrSges(3,3) * t80 + t118 * t84 + t119 * (pkin(9) * t25 + t84) + t103 * t16 - t104 * t25 + t125 * t15) * g(2) + (-m(3) * t73 + t94 * mrSges(2,1) + t24 * mrSges(3,1) + t91 * mrSges(6,1) + t93 * mrSges(7,1) + t96 * mrSges(2,2) + t90 * mrSges(6,2) + t92 * mrSges(7,2) - mrSges(3,3) * t82 + t118 * t67 + t119 * (-t23 * pkin(9) + t67) - t136 * t12 + t104 * t23 - t125 * t11) * g(1) (-t126 * (pkin(2) * t81 + pkin(9) * t88) + (((-m(6) * pkin(4) - mrSges(5,1)) * t48 - m(7) * t28 - t137 + t138) * t53 + t143 * t95) * t49) * g(3) + (t120 * t24 + t135 * t23) * g(2) + (t120 * t26 + t135 * t25) * g(1) (t125 * t22 - t136 * t21) * g(3) + (-t11 * t136 + t12 * t125) * g(2) + (t125 * t16 - t136 * t15) * g(1) (m(5) - t118) * (-g(1) * t15 - g(2) * t11 - g(3) * t21) (-(-t22 * t42 + t41 * t81) * mrSges(6,2) - t97 + t109 * (-t22 * t41 - t42 * t81)) * g(3) + (-(-t12 * t42 - t91) * mrSges(6,2) - t101 + t109 * (-t12 * t41 + t90)) * g(2) + (t8 * mrSges(6,2) + t109 * t7 - t100) * g(1), -g(1) * t100 - g(2) * t101 - g(3) * t97];
taug  = t1(:);
