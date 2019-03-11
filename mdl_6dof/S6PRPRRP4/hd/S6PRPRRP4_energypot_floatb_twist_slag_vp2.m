% Calculate potential energy for
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:08
% EndTime: 2019-03-08 20:09:09
% DurationCPUTime: 0.75s
% Computational Cost: add. (391->111), mult. (605->127), div. (0->0), fcn. (697->12), ass. (0->49)
t77 = -m(1) - m(2);
t76 = -m(6) - m(7);
t75 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t74 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t73 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t72 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t71 = -mrSges(5,3) - t73;
t41 = sin(pkin(11));
t70 = pkin(3) * t41;
t42 = sin(pkin(10));
t43 = sin(pkin(6));
t69 = t42 * t43;
t45 = cos(pkin(10));
t68 = t43 * t45;
t49 = sin(qJ(2));
t67 = t43 * t49;
t51 = cos(qJ(2));
t66 = t43 * t51;
t46 = cos(pkin(6));
t65 = t46 * t49;
t64 = t46 * t51;
t63 = t42 * pkin(1) + r_base(2);
t62 = t41 * t69;
t61 = qJ(1) + r_base(3);
t60 = t45 * pkin(1) + pkin(7) * t69 + r_base(1);
t59 = t46 * pkin(7) + t61;
t58 = -pkin(7) * t68 + t63;
t24 = t42 * t64 + t45 * t49;
t25 = -t42 * t65 + t45 * t51;
t44 = cos(pkin(11));
t34 = pkin(3) * t44 + pkin(2);
t47 = -pkin(8) - qJ(3);
t57 = pkin(3) * t62 - t24 * t47 + t25 * t34 + t60;
t56 = t34 * t67 + t46 * t70 + t47 * t66 + t59;
t22 = t42 * t49 - t45 * t64;
t23 = t42 * t51 + t45 * t65;
t53 = t23 * t34 - t22 * t47 + (-pkin(7) - t70) * t68 + t63;
t50 = cos(qJ(5));
t48 = sin(qJ(5));
t40 = pkin(11) + qJ(4);
t36 = cos(t40);
t35 = sin(t40);
t17 = t35 * t46 + t36 * t67;
t16 = t35 * t67 - t46 * t36;
t10 = t25 * t36 + t35 * t69;
t9 = t25 * t35 - t36 * t69;
t8 = t23 * t36 - t35 * t68;
t7 = t23 * t35 + t36 * t68;
t1 = (-m(1) * r_base(3) - m(2) * t61 - m(5) * t56 - t17 * mrSges(5,1) + mrSges(5,3) * t66 - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4)) * t59 + t76 * (t17 * pkin(4) + pkin(9) * t16 + t56) + (-t41 * mrSges(4,1) - t44 * mrSges(4,2) - mrSges(3,3)) * t46 + (t73 * t51 + (-m(4) * pkin(2) - t44 * mrSges(4,1) + t41 * mrSges(4,2) - mrSges(3,1)) * t49) * t43 + t74 * (t17 * t50 - t48 * t66) + t72 * (t17 * t48 + t50 * t66) + t75 * t16) * g(3) + (-m(4) * (pkin(2) * t23 + t58) - m(3) * t58 - m(5) * t53 + mrSges(3,3) * t68 - (-t23 * t41 - t44 * t68) * mrSges(4,2) - (t23 * t44 - t41 * t68) * mrSges(4,1) - mrSges(1,2) - t8 * mrSges(5,1) - t23 * mrSges(3,1) - t42 * mrSges(2,1) - t45 * mrSges(2,2) + t77 * r_base(2) + t76 * (t8 * pkin(4) + pkin(9) * t7 + t53) + t74 * (t22 * t48 + t50 * t8) + t72 * (-t22 * t50 + t48 * t8) + t75 * t7 + t71 * t22) * g(2) + (-(t25 * t44 + t62) * mrSges(4,1) - m(4) * (pkin(2) * t25 + t60) - m(3) * t60 - m(5) * t57 - mrSges(3,3) * t69 - (-t25 * t41 + t44 * t69) * mrSges(4,2) - mrSges(1,1) - t10 * mrSges(5,1) - t25 * mrSges(3,1) + t42 * mrSges(2,2) - t45 * mrSges(2,1) + t77 * r_base(1) + t76 * (t10 * pkin(4) + pkin(9) * t9 + t57) + t74 * (t10 * t50 + t24 * t48) + t72 * (t10 * t48 - t24 * t50) + t75 * t9 + t71 * t24) * g(1);
U  = t1;
