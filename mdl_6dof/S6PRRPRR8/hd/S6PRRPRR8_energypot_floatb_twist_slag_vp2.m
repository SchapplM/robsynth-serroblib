% Calculate potential energy for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:59
% EndTime: 2019-03-08 22:35:59
% DurationCPUTime: 0.80s
% Computational Cost: add. (498->109), mult. (1203->135), div. (0->0), fcn. (1483->14), ass. (0->57)
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t49 = cos(pkin(6));
t46 = sin(pkin(6));
t56 = cos(qJ(2));
t79 = t46 * t56;
t28 = -t45 * t79 + t48 * t49;
t44 = sin(pkin(12));
t47 = cos(pkin(12));
t53 = sin(qJ(2));
t75 = t49 * t56;
t31 = -t44 * t75 - t47 * t53;
t81 = t46 * t48;
t22 = -t31 * t45 + t44 * t81;
t93 = -m(1) - m(2);
t92 = -m(6) - m(7);
t91 = mrSges(4,2) - mrSges(5,3);
t90 = -mrSges(4,3) - mrSges(5,1);
t89 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t50 = sin(qJ(6));
t54 = cos(qJ(6));
t88 = -t50 * mrSges(7,1) - t54 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t87 = -m(7) * pkin(5) - mrSges(7,1) * t54 + mrSges(7,2) * t50 - mrSges(6,1);
t86 = cos(qJ(3));
t84 = t44 * t46;
t52 = sin(qJ(3));
t83 = t45 * t52;
t82 = t46 * t47;
t80 = t46 * t53;
t78 = t48 * t52;
t76 = t49 * t53;
t72 = qJ(1) + r_base(3);
t71 = t45 * t86;
t70 = t48 * t86;
t69 = pkin(1) * t47 + pkin(8) * t84 + r_base(1);
t68 = pkin(8) * t49 + t72;
t67 = t46 * t71;
t29 = -t44 * t53 + t47 * t75;
t21 = -t29 * t45 - t47 * t81;
t66 = pkin(1) * t44 - pkin(8) * t82 + r_base(2);
t32 = -t44 * t76 + t47 * t56;
t65 = t32 * pkin(2) + pkin(9) * t22 + t69;
t64 = pkin(2) * t80 + pkin(9) * t28 + t68;
t11 = -t31 * t70 + t32 * t52 - t44 * t67;
t12 = t32 * t86 + (t31 * t48 + t45 * t84) * t52;
t63 = pkin(3) * t12 + qJ(4) * t11 + t65;
t19 = -t49 * t71 + t52 * t80 - t70 * t79;
t20 = t49 * t83 + (t53 * t86 + t56 * t78) * t46;
t62 = pkin(3) * t20 + qJ(4) * t19 + t64;
t30 = t44 * t56 + t47 * t76;
t61 = pkin(2) * t30 + pkin(9) * t21 + t66;
t10 = t29 * t78 + t30 * t86 - t82 * t83;
t9 = -t29 * t70 + t30 * t52 + t47 * t67;
t58 = pkin(3) * t10 + qJ(4) * t9 + t61;
t55 = cos(qJ(5));
t51 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t72 - mrSges(2,3) - m(3) * t68 - t49 * mrSges(3,3) - (t53 * mrSges(3,1) + t56 * mrSges(3,2)) * t46 - m(4) * t64 - m(5) * t62 + t92 * (t28 * pkin(4) + t20 * pkin(10) + t62) + t90 * t28 + t91 * t19 + t89 * (-t19 * t55 + t28 * t51) + t87 * (t19 * t51 + t28 * t55) + t88 * t20) * g(3) + (-m(3) * t66 - m(4) * t61 - m(5) * t58 - t44 * mrSges(2,1) - t30 * mrSges(3,1) - t47 * mrSges(2,2) - t29 * mrSges(3,2) + mrSges(3,3) * t82 - mrSges(1,2) + t93 * r_base(2) + t91 * t9 + t92 * (t21 * pkin(4) + pkin(10) * t10 + t58) + t90 * t21 + t89 * (t21 * t51 - t9 * t55) + t87 * (t21 * t55 + t51 * t9) + t88 * t10) * g(2) + (-m(3) * t69 - m(4) * t65 - m(5) * t63 - t47 * mrSges(2,1) - t32 * mrSges(3,1) + t44 * mrSges(2,2) - t31 * mrSges(3,2) - mrSges(3,3) * t84 - mrSges(1,1) + t93 * r_base(1) + t92 * (t22 * pkin(4) + pkin(10) * t12 + t63) + t89 * (-t11 * t55 + t22 * t51) + t90 * t22 + t91 * t11 + t87 * (t11 * t51 + t22 * t55) + t88 * t12) * g(1);
U  = t1;
