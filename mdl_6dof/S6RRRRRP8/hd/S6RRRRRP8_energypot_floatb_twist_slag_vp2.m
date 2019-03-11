% Calculate potential energy for
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:09
% EndTime: 2019-03-10 01:47:10
% DurationCPUTime: 0.74s
% Computational Cost: add. (391->111), mult. (605->125), div. (0->0), fcn. (697->12), ass. (0->51)
t79 = -m(1) - m(2);
t78 = -m(6) - m(7);
t77 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t76 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3);
t75 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t74 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t73 = -mrSges(5,3) - t76;
t44 = sin(qJ(3));
t72 = pkin(3) * t44;
t41 = sin(pkin(6));
t45 = sin(qJ(2));
t71 = t41 * t45;
t46 = sin(qJ(1));
t70 = t41 * t46;
t49 = cos(qJ(2));
t69 = t41 * t49;
t50 = cos(qJ(1));
t68 = t41 * t50;
t67 = t45 * t46;
t66 = t45 * t50;
t65 = t46 * t49;
t64 = t49 * t50;
t63 = pkin(7) + r_base(3);
t62 = t46 * pkin(1) + r_base(2);
t61 = t44 * t70;
t42 = cos(pkin(6));
t60 = t42 * pkin(8) + t63;
t59 = t50 * pkin(1) + pkin(8) * t70 + r_base(1);
t58 = -pkin(8) * t68 + t62;
t48 = cos(qJ(3));
t34 = pkin(3) * t48 + pkin(2);
t51 = -pkin(10) - pkin(9);
t57 = t34 * t71 + t42 * t72 + t51 * t69 + t60;
t24 = t42 * t65 + t66;
t25 = -t42 * t67 + t64;
t56 = pkin(3) * t61 - t24 * t51 + t25 * t34 + t59;
t22 = -t42 * t64 + t67;
t23 = t42 * t66 + t65;
t53 = t23 * t34 - t22 * t51 + (-pkin(8) - t72) * t68 + t62;
t47 = cos(qJ(5));
t43 = sin(qJ(5));
t40 = qJ(3) + qJ(4);
t36 = cos(t40);
t35 = sin(t40);
t17 = t35 * t42 + t36 * t71;
t16 = t35 * t71 - t42 * t36;
t12 = t25 * t36 + t35 * t70;
t11 = t25 * t35 - t36 * t70;
t10 = t23 * t36 - t35 * t68;
t9 = t23 * t35 + t36 * t68;
t1 = (-m(1) * r_base(3) - m(2) * t63 - m(5) * t57 - t17 * mrSges(5,1) + mrSges(5,3) * t69 - mrSges(1,3) - mrSges(2,3) + t78 * (t17 * pkin(4) + pkin(11) * t16 + t57) + t75 * (t17 * t47 - t43 * t69) + t74 * (t17 * t43 + t47 * t69) + (-m(3) - m(4)) * t60 + (-t44 * mrSges(4,1) - t48 * mrSges(4,2) - mrSges(3,3)) * t42 + (t76 * t49 + (-m(4) * pkin(2) - t48 * mrSges(4,1) + t44 * mrSges(4,2) - mrSges(3,1)) * t45) * t41 + t77 * t16) * g(3) + (-m(4) * (pkin(2) * t23 + t58) - m(5) * t53 - m(3) * t58 + mrSges(3,3) * t68 - (-t23 * t44 - t48 * t68) * mrSges(4,2) - (t23 * t48 - t44 * t68) * mrSges(4,1) - mrSges(1,2) - t10 * mrSges(5,1) - t23 * mrSges(3,1) - t46 * mrSges(2,1) - t50 * mrSges(2,2) + t79 * r_base(2) + t78 * (t10 * pkin(4) + pkin(11) * t9 + t53) + t75 * (t10 * t47 + t22 * t43) + t74 * (t10 * t43 - t22 * t47) + t77 * t9 + t73 * t22) * g(2) + (-(t25 * t48 + t61) * mrSges(4,1) - m(4) * (pkin(2) * t25 + t59) - m(5) * t56 - m(3) * t59 - mrSges(3,3) * t70 - (-t25 * t44 + t48 * t70) * mrSges(4,2) - mrSges(1,1) - t12 * mrSges(5,1) - t25 * mrSges(3,1) + t46 * mrSges(2,2) - t50 * mrSges(2,1) + t79 * r_base(1) + t78 * (t12 * pkin(4) + pkin(11) * t11 + t56) + t75 * (t12 * t47 + t24 * t43) + t74 * (t12 * t43 - t24 * t47) + t73 * t24 + t77 * t11) * g(1);
U  = t1;
