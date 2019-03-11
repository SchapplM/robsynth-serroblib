% Calculate potential energy for
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:30
% EndTime: 2019-03-09 04:21:31
% DurationCPUTime: 0.80s
% Computational Cost: add. (498->109), mult. (1203->134), div. (0->0), fcn. (1483->14), ass. (0->58)
t44 = sin(pkin(12));
t49 = cos(pkin(6));
t56 = cos(qJ(1));
t47 = cos(pkin(12));
t53 = sin(qJ(1));
t76 = t53 * t47;
t31 = -t44 * t56 - t49 * t76;
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t46 = sin(pkin(6));
t82 = t46 * t53;
t22 = -t31 * t45 + t48 * t82;
t83 = t46 * t47;
t28 = -t45 * t83 + t48 * t49;
t94 = -m(1) - m(2);
t93 = -m(6) - m(7);
t92 = mrSges(4,2) - mrSges(5,3);
t91 = -mrSges(4,3) - mrSges(5,1);
t90 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t50 = sin(qJ(6));
t54 = cos(qJ(6));
t89 = -t50 * mrSges(7,1) - t54 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t88 = -m(7) * pkin(5) - t54 * mrSges(7,1) + mrSges(7,2) * t50 - mrSges(6,1);
t87 = cos(qJ(3));
t85 = t44 * t46;
t52 = sin(qJ(3));
t84 = t45 * t52;
t81 = t46 * t56;
t79 = t48 * t52;
t78 = t49 * t56;
t77 = t53 * t44;
t75 = qJ(2) * t46;
t74 = pkin(8) + r_base(3);
t71 = t45 * t87;
t70 = t48 * t87;
t69 = t49 * qJ(2) + t74;
t68 = t56 * pkin(1) + t53 * t75 + r_base(1);
t67 = t46 * t71;
t29 = t47 * t78 - t77;
t21 = -t29 * t45 - t48 * t81;
t66 = pkin(1) * t53 - t56 * t75 + r_base(2);
t32 = t47 * t56 - t49 * t77;
t65 = t32 * pkin(2) + pkin(9) * t22 + t68;
t64 = pkin(2) * t85 + pkin(9) * t28 + t69;
t13 = -t31 * t70 + t32 * t52 - t53 * t67;
t14 = t32 * t87 + (t31 * t48 + t45 * t82) * t52;
t63 = pkin(3) * t14 + qJ(4) * t13 + t65;
t17 = -t49 * t71 + t52 * t85 - t70 * t83;
t18 = t49 * t84 + (t44 * t87 + t47 * t79) * t46;
t62 = pkin(3) * t18 + qJ(4) * t17 + t64;
t30 = t44 * t78 + t76;
t61 = pkin(2) * t30 + pkin(9) * t21 + t66;
t11 = -t29 * t70 + t30 * t52 + t56 * t67;
t12 = t29 * t79 + t30 * t87 - t81 * t84;
t58 = pkin(3) * t12 + t11 * qJ(4) + t61;
t55 = cos(qJ(5));
t51 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t74 - mrSges(2,3) - m(3) * t69 - t49 * mrSges(3,3) - (t44 * mrSges(3,1) + t47 * mrSges(3,2)) * t46 - m(4) * t64 - m(5) * t62 + t93 * (t28 * pkin(4) + pkin(10) * t18 + t62) + t90 * (-t17 * t55 + t28 * t51) + t91 * t28 + t92 * t17 + t88 * (t17 * t51 + t28 * t55) + t89 * t18) * g(3) + (-m(3) * t66 - m(4) * t61 - m(5) * t58 - t53 * mrSges(2,1) - t30 * mrSges(3,1) - t56 * mrSges(2,2) - t29 * mrSges(3,2) + mrSges(3,3) * t81 - mrSges(1,2) + t94 * r_base(2) + t93 * (t21 * pkin(4) + t12 * pkin(10) + t58) + t91 * t21 + t92 * t11 + t90 * (-t11 * t55 + t21 * t51) + t88 * (t11 * t51 + t21 * t55) + t89 * t12) * g(2) + (-m(3) * t68 - m(4) * t65 - m(5) * t63 - t56 * mrSges(2,1) - t32 * mrSges(3,1) + t53 * mrSges(2,2) - t31 * mrSges(3,2) - mrSges(3,3) * t82 - mrSges(1,1) + t94 * r_base(1) + t93 * (t22 * pkin(4) + pkin(10) * t14 + t63) + t90 * (-t13 * t55 + t22 * t51) + t91 * t22 + t92 * t13 + t88 * (t13 * t51 + t22 * t55) + t89 * t14) * g(1);
U  = t1;
