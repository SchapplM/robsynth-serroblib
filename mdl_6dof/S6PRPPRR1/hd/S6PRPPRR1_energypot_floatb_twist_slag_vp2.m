% Calculate potential energy for
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:54
% EndTime: 2019-03-08 19:13:55
% DurationCPUTime: 0.83s
% Computational Cost: add. (404->97), mult. (738->111), div. (0->0), fcn. (880->14), ass. (0->53)
t52 = cos(qJ(2));
t84 = t52 * mrSges(3,2);
t83 = -m(4) - m(5);
t82 = -m(6) - m(7);
t41 = sin(pkin(11));
t45 = cos(pkin(11));
t50 = sin(qJ(2));
t81 = t50 * t41 - t45 * t52;
t80 = -m(3) - m(1) - m(2);
t79 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t44 = cos(pkin(12));
t78 = -t44 * mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t47 = cos(pkin(6));
t67 = t47 * t50;
t77 = t67 * mrSges(3,1) + t47 * t84 + mrSges(2,2);
t76 = m(3) * pkin(7) - t78;
t75 = -m(3) * pkin(1) - t52 * mrSges(3,1) + t50 * mrSges(3,2) - mrSges(2,1);
t40 = sin(pkin(12));
t74 = -m(5) * pkin(3) - t44 * mrSges(5,1) + t40 * mrSges(5,2) - mrSges(4,1);
t49 = sin(qJ(6));
t51 = cos(qJ(6));
t73 = -m(7) * pkin(5) - t51 * mrSges(7,1) + t49 * mrSges(7,2) - mrSges(6,1);
t72 = m(5) * qJ(4) + t49 * mrSges(7,1) + t51 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t71 = t40 * t47;
t42 = sin(pkin(10));
t43 = sin(pkin(6));
t70 = t42 * t43;
t46 = cos(pkin(10));
t69 = t43 * t46;
t64 = t40 * t70;
t63 = t40 * t69;
t62 = qJ(1) + r_base(3);
t22 = pkin(2) * t67 + (-pkin(7) - qJ(3)) * t43;
t34 = pkin(2) * t52 + pkin(1);
t61 = t46 * t22 + t42 * t34 + r_base(2);
t60 = t47 * pkin(7) + t62;
t59 = t41 * t52 + t50 * t45;
t58 = -t22 * t42 + t46 * t34 + r_base(1);
t57 = t43 * t50 * pkin(2) + t47 * qJ(3) + t60;
t56 = t81 * t47;
t48 = -pkin(8) - qJ(4);
t39 = pkin(12) + qJ(5);
t36 = cos(t39);
t35 = sin(t39);
t33 = pkin(4) * t44 + pkin(3);
t21 = t59 * t47;
t20 = t59 * t43;
t19 = t81 * t43;
t12 = -t21 * t42 - t46 * t81;
t11 = t42 * t56 - t46 * t59;
t10 = t21 * t46 - t42 * t81;
t9 = -t42 * t59 - t46 * t56;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t62 - mrSges(2,3) - m(3) * t60 - (t50 * mrSges(3,1) + t84) * t43 - t71 * mrSges(5,1) + t83 * t57 + t82 * (pkin(4) * t71 - t19 * t48 + t20 * t33 + t57) + t74 * t20 + t78 * t47 + t79 * (t20 * t35 - t47 * t36) + t73 * (t20 * t36 + t35 * t47) - t72 * t19) * g(3) + (t63 * mrSges(5,1) - mrSges(1,2) - t77 * t46 + t75 * t42 + t80 * r_base(2) + t83 * t61 + t74 * t10 + t76 * t69 + t82 * (-pkin(4) * t63 + t10 * t33 + t9 * t48 + t61) + t79 * (t10 * t35 + t36 * t69) + t73 * (t10 * t36 - t35 * t69) + t72 * t9) * g(2) + (-t64 * mrSges(5,1) - mrSges(1,1) + t77 * t42 + t75 * t46 + t80 * r_base(1) + t83 * t58 + t74 * t12 - t76 * t70 + t82 * (pkin(4) * t64 + t11 * t48 + t12 * t33 + t58) + t79 * (t12 * t35 - t36 * t70) + t73 * (t12 * t36 + t35 * t70) + t72 * t11) * g(1);
U  = t1;
