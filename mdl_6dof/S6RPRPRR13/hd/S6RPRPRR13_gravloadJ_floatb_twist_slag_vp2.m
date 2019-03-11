% Calculate Gravitation load on the joints for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:44
% EndTime: 2019-03-09 04:21:47
% DurationCPUTime: 0.99s
% Computational Cost: add. (909->111), mult. (2450->163), div. (0->0), fcn. (3073->14), ass. (0->65)
t111 = m(6) + m(7);
t108 = m(5) + t111;
t75 = -t108 * qJ(4) + mrSges(4,2) - mrSges(5,3);
t55 = sin(qJ(6));
t59 = cos(qJ(6));
t76 = -pkin(10) * t111 - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t115 = -t55 * mrSges(7,1) - t59 * mrSges(7,2) + t76;
t110 = m(7) * pkin(5) + t59 * mrSges(7,1) - t55 * mrSges(7,2) + mrSges(6,1);
t114 = m(7) * pkin(11) - mrSges(6,2) + mrSges(7,3);
t51 = sin(pkin(7));
t58 = sin(qJ(1));
t50 = sin(pkin(12));
t54 = cos(pkin(6));
t61 = cos(qJ(1));
t53 = cos(pkin(12));
t95 = t58 * t53;
t77 = t50 * t61 + t54 * t95;
t52 = sin(pkin(6));
t91 = cos(pkin(7));
t86 = t52 * t91;
t30 = t77 * t51 + t58 * t86;
t96 = t58 * t50;
t97 = t54 * t61;
t37 = -t53 * t97 + t96;
t29 = t37 * t51 - t61 * t86;
t113 = pkin(3) * t108 - t115;
t38 = t50 * t97 + t95;
t57 = sin(qJ(3));
t102 = cos(qJ(3));
t79 = t91 * t102;
t88 = t51 * t102;
t82 = t52 * t88;
t17 = t37 * t79 + t38 * t57 + t61 * t82;
t56 = sin(qJ(5));
t60 = cos(qJ(5));
t112 = t17 * t60 - t29 * t56;
t6 = t17 * t56 + t29 * t60;
t104 = -t110 * t56 + t114 * t60 + t75;
t100 = t52 * t53;
t99 = t52 * t57;
t98 = t52 * t58;
t94 = -mrSges(4,3) - mrSges(5,1);
t92 = qJ(2) * t52;
t93 = t61 * pkin(1) + t58 * t92;
t87 = -pkin(1) * t58 + t61 * t92;
t85 = t57 * t91;
t69 = -t38 * pkin(2) - t29 * pkin(9) + t87;
t68 = t77 * t91;
t18 = -t61 * t51 * t99 + t38 * t102 - t37 * t85;
t67 = -pkin(3) * t18 + t69;
t39 = t53 * t61 - t54 * t96;
t65 = t39 * pkin(2) + t30 * pkin(9) + t93;
t22 = t39 * t102 + (t51 * t98 - t68) * t57;
t64 = t22 * pkin(3) + t65;
t63 = t30 * pkin(4) + t64;
t36 = -t51 * t100 + t54 * t91;
t25 = t54 * t51 * t57 + (t102 * t50 + t53 * t85) * t52;
t24 = -t79 * t100 + t50 * t99 - t54 * t88;
t21 = t102 * t68 + t39 * t57 - t58 * t82;
t16 = t24 * t56 + t36 * t60;
t8 = t21 * t56 + t30 * t60;
t7 = -t21 * t60 + t30 * t56;
t2 = t22 * t55 + t59 * t8;
t1 = t22 * t59 - t55 * t8;
t3 = [(-t61 * mrSges(2,1) + t58 * mrSges(2,2) - m(3) * t93 - t39 * mrSges(3,1) + t77 * mrSges(3,2) - mrSges(3,3) * t98 - m(4) * t65 - m(5) * t64 - m(6) * t63 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t63) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t114 * t7 + t94 * t30 + t75 * t21 + t76 * t22) * g(2) + (t58 * mrSges(2,1) - m(3) * t87 + t38 * mrSges(3,1) - t37 * mrSges(3,2) - m(4) * t69 - m(5) * t67 + (-t52 * mrSges(3,3) + mrSges(2,2)) * t61 - t114 * t112 - t94 * t29 - t75 * t17 + t110 * t6 - t115 * t18 + t111 * (pkin(4) * t29 - t67)) * g(1) (-t54 * g(3) + (-t58 * g(1) + t61 * g(2)) * t52) * (m(3) + m(4) + t108) (t104 * t25 + t113 * t24) * g(3) + (t104 * t18 + t113 * t17) * g(2) + (t104 * t22 + t113 * t21) * g(1), t108 * (-g(1) * t21 - g(2) * t17 - g(3) * t24) (-t114 * t16 - t110 * (t24 * t60 - t36 * t56)) * g(3) + (-t110 * t112 - t114 * t6) * g(2) + (t110 * t7 - t114 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t18 * t59 - t55 * t6) * mrSges(7,1) + (-t18 * t55 - t59 * t6) * mrSges(7,2)) - g(3) * ((-t16 * t55 + t25 * t59) * mrSges(7,1) + (-t16 * t59 - t25 * t55) * mrSges(7,2))];
taug  = t3(:);
