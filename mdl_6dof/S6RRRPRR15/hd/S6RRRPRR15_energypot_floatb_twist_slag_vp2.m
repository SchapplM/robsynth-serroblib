% Calculate potential energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR15_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:25
% EndTime: 2019-03-09 20:24:26
% DurationCPUTime: 0.82s
% Computational Cost: add. (498->109), mult. (1203->132), div. (0->0), fcn. (1483->14), ass. (0->58)
t47 = cos(pkin(6));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t76 = t52 * t55;
t51 = sin(qJ(2));
t56 = cos(qJ(1));
t78 = t51 * t56;
t31 = -t47 * t76 - t78;
t44 = sin(pkin(7));
t46 = cos(pkin(7));
t45 = sin(pkin(6));
t83 = t45 * t52;
t22 = -t31 * t44 + t46 * t83;
t82 = t45 * t55;
t28 = -t44 * t82 + t46 * t47;
t94 = -m(1) - m(2);
t93 = -m(6) - m(7);
t92 = mrSges(4,2) - mrSges(5,3);
t91 = -mrSges(4,3) - mrSges(5,1);
t90 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t48 = sin(qJ(6));
t53 = cos(qJ(6));
t89 = -t48 * mrSges(7,1) - t53 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t88 = -m(7) * pkin(5) - t53 * mrSges(7,1) + mrSges(7,2) * t48 - mrSges(6,1);
t87 = cos(qJ(3));
t50 = sin(qJ(3));
t85 = t44 * t50;
t84 = t45 * t51;
t81 = t45 * t56;
t79 = t46 * t50;
t77 = t52 * t51;
t75 = t55 * t56;
t74 = pkin(8) + r_base(3);
t71 = t44 * t87;
t70 = t46 * t87;
t69 = pkin(9) * t47 + t74;
t68 = pkin(1) * t56 + pkin(9) * t83 + r_base(1);
t67 = t45 * t71;
t29 = t47 * t75 - t77;
t21 = -t29 * t44 - t46 * t81;
t66 = pkin(1) * t52 - pkin(9) * t81 + r_base(2);
t32 = -t47 * t77 + t75;
t65 = t32 * pkin(2) + pkin(10) * t22 + t68;
t64 = pkin(2) * t84 + pkin(10) * t28 + t69;
t13 = -t31 * t70 + t32 * t50 - t52 * t67;
t14 = t32 * t87 + (t31 * t46 + t44 * t83) * t50;
t63 = pkin(3) * t14 + qJ(4) * t13 + t65;
t17 = -t47 * t71 + t50 * t84 - t70 * t82;
t18 = t47 * t85 + (t51 * t87 + t55 * t79) * t45;
t62 = pkin(3) * t18 + qJ(4) * t17 + t64;
t30 = t47 * t78 + t76;
t61 = pkin(2) * t30 + pkin(10) * t21 + t66;
t11 = -t29 * t70 + t30 * t50 + t56 * t67;
t12 = t29 * t79 + t30 * t87 - t81 * t85;
t58 = pkin(3) * t12 + t11 * qJ(4) + t61;
t54 = cos(qJ(5));
t49 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t74 - mrSges(2,3) - m(3) * t69 - t47 * mrSges(3,3) - (t51 * mrSges(3,1) + t55 * mrSges(3,2)) * t45 - m(4) * t64 - m(5) * t62 + t93 * (t28 * pkin(4) + pkin(11) * t18 + t62) + t90 * (-t17 * t54 + t28 * t49) + t91 * t28 + t92 * t17 + t88 * (t17 * t49 + t28 * t54) + t89 * t18) * g(3) + (-m(3) * t66 - m(4) * t61 - m(5) * t58 - t52 * mrSges(2,1) - t30 * mrSges(3,1) - t56 * mrSges(2,2) - t29 * mrSges(3,2) + mrSges(3,3) * t81 - mrSges(1,2) + t94 * r_base(2) + t93 * (t21 * pkin(4) + t12 * pkin(11) + t58) + t91 * t21 + t92 * t11 + t90 * (-t11 * t54 + t21 * t49) + t88 * (t11 * t49 + t21 * t54) + t89 * t12) * g(2) + (-m(3) * t68 - m(4) * t65 - m(5) * t63 - t56 * mrSges(2,1) - t32 * mrSges(3,1) + t52 * mrSges(2,2) - t31 * mrSges(3,2) - mrSges(3,3) * t83 - mrSges(1,1) + t94 * r_base(1) + t93 * (t22 * pkin(4) + pkin(11) * t14 + t63) + t90 * (-t13 * t54 + t22 * t49) + t91 * t22 + t92 * t13 + t88 * (t13 * t49 + t22 * t54) + t89 * t14) * g(1);
U  = t1;
