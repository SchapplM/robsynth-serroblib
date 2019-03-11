% Calculate potential energy for
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energypot_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:17
% EndTime: 2019-03-09 07:50:18
% DurationCPUTime: 0.93s
% Computational Cost: add. (972->128), mult. (2579->169), div. (0->0), fcn. (3324->18), ass. (0->68)
t58 = sin(pkin(7));
t59 = sin(pkin(6));
t60 = cos(pkin(14));
t62 = cos(pkin(7));
t63 = cos(pkin(6));
t44 = -t58 * t59 * t60 + t62 * t63;
t56 = sin(pkin(14));
t72 = cos(qJ(1));
t68 = sin(qJ(1));
t93 = t68 * t60;
t47 = -t56 * t72 - t63 * t93;
t99 = t59 * t68;
t39 = -t47 * t58 + t62 * t99;
t100 = t58 * t63;
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t97 = t60 * t62;
t36 = t71 * t100 + (-t56 * t67 + t71 * t97) * t59;
t57 = sin(pkin(8));
t61 = cos(pkin(8));
t23 = -t36 * t57 + t44 * t61;
t94 = t68 * t56;
t48 = t60 * t72 - t63 * t94;
t83 = t47 * t62 + t58 * t99;
t28 = -t48 * t67 + t71 * t83;
t19 = -t28 * t57 + t39 * t61;
t95 = t63 * t72;
t46 = t56 * t95 + t93;
t45 = t60 * t95 - t94;
t98 = t59 * t72;
t84 = t45 * t62 - t58 * t98;
t26 = -t46 * t67 + t71 * t84;
t38 = -t45 * t58 - t62 * t98;
t18 = -t26 * t57 + t38 * t61;
t113 = -m(1) - m(2);
t112 = -m(6) - m(7);
t111 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t64 = sin(qJ(6));
t69 = cos(qJ(6));
t110 = -t64 * mrSges(7,1) - t69 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t109 = -m(7) * pkin(5) - t69 * mrSges(7,1) + t64 * mrSges(7,2) - mrSges(6,1);
t108 = cos(qJ(4));
t92 = qJ(2) * t59;
t91 = pkin(9) + r_base(3);
t88 = t57 * t108;
t87 = t61 * t108;
t86 = t63 * qJ(2) + t91;
t85 = t72 * pkin(1) + t68 * t92 + r_base(1);
t82 = t68 * pkin(1) - t72 * t92 + r_base(2);
t81 = t48 * pkin(2) + t39 * pkin(10) + t85;
t80 = t59 * t56 * pkin(2) + t44 * pkin(10) + t86;
t29 = t48 * t71 + t67 * t83;
t79 = t29 * pkin(3) + t19 * pkin(11) + t81;
t37 = t67 * t100 + (t56 * t71 + t67 * t97) * t59;
t78 = t37 * pkin(3) + t23 * pkin(11) + t80;
t77 = t46 * pkin(2) + pkin(10) * t38 + t82;
t27 = t46 * t71 + t67 * t84;
t74 = t27 * pkin(3) + t18 * pkin(11) + t77;
t70 = cos(qJ(5));
t66 = sin(qJ(4));
t65 = sin(qJ(5));
t15 = t37 * t108 + (t36 * t61 + t44 * t57) * t66;
t14 = -t36 * t87 + t37 * t66 - t44 * t88;
t12 = t29 * t108 + (t28 * t61 + t39 * t57) * t66;
t11 = -t28 * t87 + t29 * t66 - t39 * t88;
t10 = t27 * t108 + (t26 * t61 + t38 * t57) * t66;
t9 = -t26 * t87 + t27 * t66 - t38 * t88;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t91 - mrSges(2,3) - m(3) * t86 - t63 * mrSges(3,3) - (t56 * mrSges(3,1) + t60 * mrSges(3,2)) * t59 - m(4) * t80 - t37 * mrSges(4,1) - t36 * mrSges(4,2) - t44 * mrSges(4,3) - m(5) * t78 - t15 * mrSges(5,1) - t23 * mrSges(5,3) + t112 * (t15 * pkin(4) + pkin(12) * t14 + t78) + t111 * (t15 * t65 - t23 * t70) + t109 * (t15 * t70 + t23 * t65) + t110 * t14) * g(3) + (-m(3) * t82 - m(4) * t77 - m(5) * t74 - t68 * mrSges(2,1) - t46 * mrSges(3,1) - t27 * mrSges(4,1) - t10 * mrSges(5,1) - t72 * mrSges(2,2) - t45 * mrSges(3,2) - t26 * mrSges(4,2) + mrSges(3,3) * t98 - t38 * mrSges(4,3) - t18 * mrSges(5,3) - mrSges(1,2) + t113 * r_base(2) + t112 * (t10 * pkin(4) + t9 * pkin(12) + t74) + t109 * (t10 * t70 + t18 * t65) + t110 * t9 + t111 * (t10 * t65 - t18 * t70)) * g(2) + (-m(3) * t85 - m(4) * t81 - m(5) * t79 - t72 * mrSges(2,1) - t48 * mrSges(3,1) - t29 * mrSges(4,1) - t12 * mrSges(5,1) + t68 * mrSges(2,2) - t47 * mrSges(3,2) - t28 * mrSges(4,2) - mrSges(3,3) * t99 - t39 * mrSges(4,3) - t19 * mrSges(5,3) - mrSges(1,1) + t113 * r_base(1) + t112 * (t12 * pkin(4) + pkin(12) * t11 + t79) + t111 * (t12 * t65 - t19 * t70) + t109 * (t12 * t70 + t19 * t65) + t110 * t11) * g(1);
U  = t1;
