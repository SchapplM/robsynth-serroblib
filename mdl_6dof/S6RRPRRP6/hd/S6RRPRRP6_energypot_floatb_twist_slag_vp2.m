% Calculate potential energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:39
% EndTime: 2019-03-09 12:07:39
% DurationCPUTime: 0.68s
% Computational Cost: add. (439->105), mult. (952->128), div. (0->0), fcn. (1177->12), ass. (0->56)
t84 = -m(6) - m(7);
t83 = mrSges(4,2) - mrSges(5,3);
t82 = -mrSges(3,3) - mrSges(4,3);
t81 = -m(3) - m(1) - m(2);
t80 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t79 = -m(3) * pkin(1) - mrSges(2,1);
t78 = m(3) * pkin(8) - t82;
t77 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t76 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t44 = sin(pkin(6));
t49 = sin(qJ(2));
t75 = t44 * t49;
t50 = sin(qJ(1));
t74 = t44 * t50;
t54 = cos(qJ(1));
t73 = t44 * t54;
t45 = cos(pkin(11));
t53 = cos(qJ(2));
t72 = t45 * t53;
t71 = t49 * t54;
t70 = t50 * t49;
t69 = t50 * t53;
t68 = t53 * t54;
t67 = pkin(7) + r_base(3);
t46 = cos(pkin(6));
t66 = t46 * pkin(8) + t67;
t31 = pkin(2) * t46 * t49 + (-pkin(8) - qJ(3)) * t44;
t40 = pkin(2) * t53 + pkin(1);
t65 = t54 * t31 + t50 * t40 + r_base(2);
t43 = sin(pkin(11));
t64 = t43 * t53 + t45 * t49;
t33 = -t43 * t49 + t72;
t63 = -t31 * t50 + t54 * t40 + r_base(1);
t62 = pkin(2) * t75 + t46 * qJ(3) + t66;
t61 = t33 * t46;
t17 = -t50 * t64 + t54 * t61;
t30 = t64 * t46;
t18 = t30 * t54 + t50 * t33;
t60 = t18 * pkin(3) - pkin(9) * t17 + t65;
t19 = -t50 * t61 - t54 * t64;
t20 = -t50 * t30 + t33 * t54;
t59 = t20 * pkin(3) - pkin(9) * t19 + t63;
t28 = t43 * t75 - t44 * t72;
t29 = t64 * t44;
t58 = t29 * pkin(3) + pkin(9) * t28 + t62;
t52 = cos(qJ(4));
t51 = cos(qJ(5));
t48 = sin(qJ(4));
t47 = sin(qJ(5));
t23 = t29 * t52 + t46 * t48;
t22 = t29 * t48 - t46 * t52;
t12 = t20 * t52 + t48 * t74;
t11 = t20 * t48 - t52 * t74;
t10 = t18 * t52 - t48 * t73;
t9 = t18 * t48 + t52 * t73;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t67 - mrSges(2,3) - m(3) * t66 - (t49 * mrSges(3,1) + t53 * mrSges(3,2)) * t44 - m(4) * t62 - t29 * mrSges(4,1) - m(5) * t58 - t23 * mrSges(5,1) + t84 * (t23 * pkin(4) + pkin(10) * t22 + t58) + t77 * (t23 * t51 + t28 * t47) + t76 * (t23 * t47 - t28 * t51) + t82 * t46 + t83 * t28 + t80 * t22) * g(3) + (-(t46 * t68 - t70) * mrSges(3,2) - (t46 * t71 + t69) * mrSges(3,1) - m(4) * t65 - m(5) * t60 - t54 * mrSges(2,2) - t18 * mrSges(4,1) - t10 * mrSges(5,1) - mrSges(1,2) + t79 * t50 + t81 * r_base(2) + t78 * t73 + t84 * (t10 * pkin(4) + pkin(10) * t9 + t60) + t77 * (t10 * t51 - t17 * t47) - t83 * t17 + t76 * (t10 * t47 + t17 * t51) + t80 * t9) * g(2) + (-(-t46 * t70 + t68) * mrSges(3,1) - (-t46 * t69 - t71) * mrSges(3,2) - m(4) * t63 - m(5) * t59 + t50 * mrSges(2,2) - t20 * mrSges(4,1) - t12 * mrSges(5,1) - mrSges(1,1) + t79 * t54 + t81 * r_base(1) - t78 * t74 + t84 * (t12 * pkin(4) + pkin(10) * t11 + t59) + t77 * (t12 * t51 - t19 * t47) + t76 * (t12 * t47 + t19 * t51) - t83 * t19 + t80 * t11) * g(1);
U  = t1;
