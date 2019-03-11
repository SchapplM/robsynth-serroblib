% Calculate Gravitation load on the joints for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:34
% EndTime: 2019-03-09 08:54:37
% DurationCPUTime: 0.98s
% Computational Cost: add. (786->116), mult. (1640->162), div. (0->0), fcn. (1998->14), ass. (0->62)
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t115 = m(7) * pkin(5) + t63 * mrSges(7,1) - t60 * mrSges(7,2) + mrSges(6,1);
t84 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t116 = m(6) + m(7);
t113 = m(5) + t116;
t91 = m(4) + t113;
t108 = cos(qJ(2));
t61 = sin(qJ(2));
t97 = sin(pkin(11));
t98 = cos(pkin(11));
t69 = t108 * t97 + t61 * t98;
t54 = pkin(12) + qJ(5);
t52 = sin(t54);
t53 = cos(t54);
t55 = sin(pkin(12));
t57 = cos(pkin(12));
t71 = m(5) * pkin(3) + t57 * mrSges(5,1) - t55 * mrSges(5,2) + mrSges(4,1);
t110 = -t115 * t53 + t84 * t52 - t71;
t62 = sin(qJ(1));
t102 = t62 * t61;
t58 = cos(pkin(6));
t64 = cos(qJ(1));
t87 = t64 * t108;
t32 = -t58 * t87 + t102;
t117 = -m(4) - m(5);
t76 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t114 = -mrSges(7,1) * t60 - mrSges(7,2) * t63 - t76;
t36 = -t108 * t98 + t61 * t97;
t100 = t64 * t61;
t88 = t62 * t108;
t34 = -t58 * t88 - t100;
t92 = t108 * pkin(2);
t111 = t91 * (t92 + pkin(1)) + m(3) * pkin(1) + mrSges(2,1);
t99 = t69 * t58;
t22 = -t64 * t36 - t62 * t99;
t17 = t62 * t36 - t64 * t99;
t56 = sin(pkin(6));
t105 = t56 * t62;
t104 = t56 * t64;
t96 = pkin(4) * t55 * t56;
t67 = t58 * t36;
t21 = t62 * t67 - t64 * t69;
t50 = pkin(4) * t57 + pkin(3);
t59 = -pkin(9) - qJ(4);
t94 = t21 * t59 + t22 * t50 + t62 * t96;
t4 = -t52 * t104 - t17 * t53;
t3 = -t53 * t104 + t17 * t52;
t48 = t56 * t92;
t81 = t32 * pkin(2);
t65 = mrSges(2,2) + t91 * (pkin(2) * t58 * t61 + (-pkin(8) - qJ(3)) * t56) + (-m(3) * pkin(8) - mrSges(5,1) * t55 - mrSges(5,2) * t57 - mrSges(3,3) - mrSges(4,3)) * t56;
t35 = -t58 * t102 + t87;
t33 = -t58 * t100 - t88;
t30 = t69 * t56;
t29 = t36 * t56;
t24 = t30 * t53 + t52 * t58;
t18 = -t62 * t69 - t64 * t67;
t8 = t52 * t105 + t22 * t53;
t7 = -t53 * t105 + t22 * t52;
t2 = -t21 * t60 + t63 * t8;
t1 = -t21 * t63 - t60 * t8;
t5 = [(-t35 * mrSges(3,1) - t34 * mrSges(3,2) - m(6) * t94 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t94) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t84 * t7 - t71 * t22 + t65 * t62 + t76 * t21 - t111 * t64) * g(2) + (-t33 * mrSges(3,1) - t32 * mrSges(3,2) + t84 * t3 + t111 * t62 - t71 * t17 + t65 * t64 + t115 * t4 + t114 * t18 + t116 * (-t17 * t50 + t18 * t59 - t64 * t96)) * g(1) (-(mrSges(3,1) * t108 - mrSges(3,2) * t61) * t56 - t116 * (-t29 * t50 - t30 * t59 + t48) + t117 * t48 + t114 * t30 - t110 * t29) * g(3) + (mrSges(3,1) * t32 - mrSges(3,2) * t33 - t117 * t81 - t116 * (t17 * t59 + t18 * t50 - t81) + t110 * t18 - t114 * t17) * g(2) + (mrSges(3,2) * t35 - t116 * (t21 * t50 - t22 * t59) + t114 * t22 + t110 * t21 + (-t91 * pkin(2) - mrSges(3,1)) * t34) * g(1) (-t58 * g(3) + (-t62 * g(1) + t64 * g(2)) * t56) * t91, t113 * (g(1) * t21 + g(2) * t18 - g(3) * t29) (t84 * t24 - t115 * (-t30 * t52 + t53 * t58)) * g(3) + (-t115 * t3 + t84 * t4) * g(2) + (t115 * t7 + t84 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t18 * t63 - t4 * t60) * mrSges(7,1) + (t18 * t60 - t4 * t63) * mrSges(7,2)) - g(3) * ((-t24 * t60 + t29 * t63) * mrSges(7,1) + (-t24 * t63 - t29 * t60) * mrSges(7,2))];
taug  = t5(:);
