% Calculate potential energy for
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:35
% EndTime: 2019-03-10 05:49:36
% DurationCPUTime: 0.93s
% Computational Cost: add. (972->128), mult. (2579->167), div. (0->0), fcn. (3324->18), ass. (0->68)
t57 = sin(pkin(7));
t58 = sin(pkin(6));
t60 = cos(pkin(7));
t61 = cos(pkin(6));
t71 = cos(qJ(2));
t44 = -t57 * t58 * t71 + t60 * t61;
t67 = sin(qJ(1));
t93 = t67 * t71;
t66 = sin(qJ(2));
t72 = cos(qJ(1));
t95 = t66 * t72;
t47 = -t61 * t93 - t95;
t99 = t58 * t67;
t39 = -t47 * t57 + t60 * t99;
t100 = t57 * t61;
t65 = sin(qJ(3));
t70 = cos(qJ(3));
t96 = t60 * t71;
t36 = t70 * t100 + (-t65 * t66 + t70 * t96) * t58;
t56 = sin(pkin(8));
t59 = cos(pkin(8));
t23 = -t36 * t56 + t44 * t59;
t92 = t71 * t72;
t94 = t67 * t66;
t48 = -t61 * t94 + t92;
t83 = t47 * t60 + t57 * t99;
t28 = -t48 * t65 + t70 * t83;
t19 = -t28 * t56 + t39 * t59;
t46 = t61 * t95 + t93;
t45 = t61 * t92 - t94;
t98 = t58 * t72;
t84 = t45 * t60 - t57 * t98;
t26 = -t46 * t65 + t70 * t84;
t38 = -t45 * t57 - t60 * t98;
t18 = -t26 * t56 + t38 * t59;
t113 = -m(1) - m(2);
t112 = -m(6) - m(7);
t111 = -m(7) * pkin(14) + mrSges(6,2) - mrSges(7,3);
t62 = sin(qJ(6));
t68 = cos(qJ(6));
t110 = -t62 * mrSges(7,1) - t68 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t109 = -m(7) * pkin(5) - mrSges(7,1) * t68 + mrSges(7,2) * t62 - mrSges(6,1);
t108 = cos(qJ(4));
t91 = pkin(9) + r_base(3);
t88 = t56 * t108;
t87 = t59 * t108;
t86 = t61 * pkin(10) + t91;
t85 = t72 * pkin(1) + pkin(10) * t99 + r_base(1);
t82 = t67 * pkin(1) - pkin(10) * t98 + r_base(2);
t81 = t48 * pkin(2) + t39 * pkin(11) + t85;
t80 = t58 * t66 * pkin(2) + t44 * pkin(11) + t86;
t29 = t48 * t70 + t65 * t83;
t79 = t29 * pkin(3) + t19 * pkin(12) + t81;
t37 = t65 * t100 + (t65 * t96 + t66 * t70) * t58;
t78 = t37 * pkin(3) + t23 * pkin(12) + t80;
t77 = t46 * pkin(2) + pkin(11) * t38 + t82;
t27 = t46 * t70 + t65 * t84;
t74 = t27 * pkin(3) + t18 * pkin(12) + t77;
t69 = cos(qJ(5));
t64 = sin(qJ(4));
t63 = sin(qJ(5));
t15 = t37 * t108 + (t36 * t59 + t44 * t56) * t64;
t14 = -t36 * t87 + t37 * t64 - t44 * t88;
t12 = t29 * t108 + (t28 * t59 + t39 * t56) * t64;
t11 = -t28 * t87 + t29 * t64 - t39 * t88;
t10 = t27 * t108 + (t26 * t59 + t38 * t56) * t64;
t9 = -t26 * t87 + t27 * t64 - t38 * t88;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t91 - mrSges(2,3) - m(3) * t86 - t61 * mrSges(3,3) - (t66 * mrSges(3,1) + t71 * mrSges(3,2)) * t58 - m(4) * t80 - t37 * mrSges(4,1) - t36 * mrSges(4,2) - t44 * mrSges(4,3) - m(5) * t78 - t15 * mrSges(5,1) - t23 * mrSges(5,3) + t112 * (t15 * pkin(4) + pkin(13) * t14 + t78) + t111 * (t15 * t63 - t23 * t69) + t109 * (t15 * t69 + t23 * t63) + t110 * t14) * g(3) + (-m(3) * t82 - m(4) * t77 - m(5) * t74 - t67 * mrSges(2,1) - t46 * mrSges(3,1) - t27 * mrSges(4,1) - t10 * mrSges(5,1) - t72 * mrSges(2,2) - t45 * mrSges(3,2) - t26 * mrSges(4,2) + mrSges(3,3) * t98 - t38 * mrSges(4,3) - t18 * mrSges(5,3) - mrSges(1,2) + t113 * r_base(2) + t112 * (t10 * pkin(4) + t9 * pkin(13) + t74) + t109 * (t10 * t69 + t18 * t63) + t110 * t9 + t111 * (t10 * t63 - t18 * t69)) * g(2) + (-m(3) * t85 - m(4) * t81 - m(5) * t79 - t72 * mrSges(2,1) - t48 * mrSges(3,1) - t29 * mrSges(4,1) - t12 * mrSges(5,1) + t67 * mrSges(2,2) - t47 * mrSges(3,2) - t28 * mrSges(4,2) - mrSges(3,3) * t99 - t39 * mrSges(4,3) - t19 * mrSges(5,3) - mrSges(1,1) + t113 * r_base(1) + t112 * (t12 * pkin(4) + pkin(13) * t11 + t79) + t111 * (t12 * t63 - t19 * t69) + t109 * (t12 * t69 + t19 * t63) + t110 * t11) * g(1);
U  = t1;
