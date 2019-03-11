% Calculate potential energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPPRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:02
% EndTime: 2019-03-08 18:39:03
% DurationCPUTime: 0.94s
% Computational Cost: add. (972->128), mult. (2579->169), div. (0->0), fcn. (3324->18), ass. (0->69)
t63 = cos(pkin(13));
t66 = cos(pkin(7));
t67 = cos(pkin(6));
t60 = sin(pkin(7));
t61 = sin(pkin(6));
t99 = t60 * t61;
t44 = -t63 * t99 + t66 * t67;
t58 = sin(pkin(12));
t100 = t58 * t67;
t57 = sin(pkin(13));
t64 = cos(pkin(12));
t47 = -t100 * t63 - t57 * t64;
t96 = t61 * t66;
t39 = -t47 * t60 + t58 * t96;
t56 = sin(pkin(14));
t62 = cos(pkin(14));
t95 = t63 * t66;
t98 = t60 * t67;
t36 = t62 * t98 + (-t56 * t57 + t62 * t95) * t61;
t59 = sin(pkin(8));
t65 = cos(pkin(8));
t29 = -t36 * t59 + t44 * t65;
t48 = -t100 * t57 + t63 * t64;
t83 = t47 * t66 + t58 * t99;
t27 = -t48 * t56 + t62 * t83;
t19 = -t27 * t59 + t39 * t65;
t94 = t64 * t67;
t46 = t57 * t94 + t58 * t63;
t45 = -t57 * t58 + t63 * t94;
t97 = t61 * t64;
t84 = t45 * t66 - t60 * t97;
t25 = -t46 * t56 + t62 * t84;
t38 = -t45 * t60 - t64 * t96;
t18 = -t25 * t59 + t38 * t65;
t114 = -m(1) - m(2);
t113 = -m(6) - m(7);
t112 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t111 = -t68 * mrSges(7,1) - t71 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t110 = -m(7) * pkin(5) - t71 * mrSges(7,1) + t68 * mrSges(7,2) - mrSges(6,1);
t109 = cos(qJ(4));
t101 = t57 * t61;
t92 = qJ(2) * t61;
t89 = qJ(1) + r_base(3);
t88 = t59 * t109;
t87 = t65 * t109;
t86 = t64 * pkin(1) + t58 * t92 + r_base(1);
t85 = t67 * qJ(2) + t89;
t82 = t58 * pkin(1) - t64 * t92 + r_base(2);
t81 = t48 * pkin(2) + t39 * qJ(3) + t86;
t80 = pkin(2) * t101 + t44 * qJ(3) + t85;
t28 = t48 * t62 + t56 * t83;
t79 = t28 * pkin(3) + t19 * pkin(9) + t81;
t37 = t62 * t101 + (t61 * t95 + t98) * t56;
t78 = t37 * pkin(3) + t29 * pkin(9) + t80;
t77 = t46 * pkin(2) + qJ(3) * t38 + t82;
t26 = t46 * t62 + t56 * t84;
t74 = t26 * pkin(3) + t18 * pkin(9) + t77;
t72 = cos(qJ(5));
t70 = sin(qJ(4));
t69 = sin(qJ(5));
t17 = t37 * t109 + (t36 * t65 + t44 * t59) * t70;
t16 = -t36 * t87 + t37 * t70 - t44 * t88;
t12 = t28 * t109 + (t27 * t65 + t39 * t59) * t70;
t11 = -t27 * t87 + t28 * t70 - t39 * t88;
t10 = t26 * t109 + (t25 * t65 + t38 * t59) * t70;
t9 = -t25 * t87 + t26 * t70 - t38 * t88;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t89 - mrSges(2,3) - m(3) * t85 - t67 * mrSges(3,3) - (t57 * mrSges(3,1) + t63 * mrSges(3,2)) * t61 - m(4) * t80 - t37 * mrSges(4,1) - t36 * mrSges(4,2) - t44 * mrSges(4,3) - m(5) * t78 - t17 * mrSges(5,1) - t29 * mrSges(5,3) + t113 * (t17 * pkin(4) + pkin(10) * t16 + t78) + t112 * (t17 * t69 - t29 * t72) + t110 * (t17 * t72 + t29 * t69) + t111 * t16) * g(3) + (-m(3) * t82 - m(4) * t77 - m(5) * t74 - t58 * mrSges(2,1) - t46 * mrSges(3,1) - t26 * mrSges(4,1) - t10 * mrSges(5,1) - t64 * mrSges(2,2) - t45 * mrSges(3,2) - t25 * mrSges(4,2) + mrSges(3,3) * t97 - t38 * mrSges(4,3) - t18 * mrSges(5,3) - mrSges(1,2) + t114 * r_base(2) + t113 * (t10 * pkin(4) + pkin(10) * t9 + t74) + t110 * (t10 * t72 + t18 * t69) + t111 * t9 + t112 * (t10 * t69 - t18 * t72)) * g(2) + (-m(3) * t86 - m(4) * t81 - m(5) * t79 - t64 * mrSges(2,1) - t48 * mrSges(3,1) - t28 * mrSges(4,1) - t12 * mrSges(5,1) - t47 * mrSges(3,2) - t27 * mrSges(4,2) - t39 * mrSges(4,3) - t19 * mrSges(5,3) - mrSges(1,1) + t114 * r_base(1) + t113 * (t12 * pkin(4) + pkin(10) * t11 + t79) + (-t61 * mrSges(3,3) + mrSges(2,2)) * t58 + t112 * (t12 * t69 - t19 * t72) + t110 * (t12 * t72 + t19 * t69) + t111 * t11) * g(1);
U  = t1;
