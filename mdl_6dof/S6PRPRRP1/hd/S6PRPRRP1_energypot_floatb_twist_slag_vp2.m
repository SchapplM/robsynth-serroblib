% Calculate potential energy for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:07
% EndTime: 2019-03-08 19:55:08
% DurationCPUTime: 0.81s
% Computational Cost: add. (409->104), mult. (867->119), div. (0->0), fcn. (1059->12), ass. (0->52)
t49 = cos(qJ(2));
t83 = t49 * mrSges(3,2);
t45 = sin(qJ(4));
t73 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t82 = t73 * t45 - mrSges(4,1);
t40 = sin(pkin(6));
t42 = cos(pkin(6));
t48 = cos(qJ(4));
t46 = sin(qJ(2));
t67 = t42 * t46;
t78 = -mrSges(3,3) - mrSges(4,3);
t81 = -t67 * mrSges(3,1) - t42 * t83 - mrSges(2,2) + (m(3) * pkin(7) + t73 * t48 - t78) * t40;
t80 = -mrSges(6,1) - mrSges(7,1);
t79 = -mrSges(6,2) - mrSges(7,2);
t38 = sin(pkin(11));
t64 = cos(pkin(11));
t27 = -t46 * t38 + t49 * t64;
t77 = -m(3) - m(2) - m(1);
t44 = sin(qJ(5));
t75 = m(7) * (pkin(5) * t44 + pkin(8)) - mrSges(4,2) + mrSges(5,3);
t72 = -m(3) * pkin(1) - t49 * mrSges(3,1) + t46 * mrSges(3,2) - mrSges(2,1);
t69 = t40 * t45;
t63 = qJ(1) + r_base(3);
t25 = pkin(2) * t67 + (-pkin(7) - qJ(3)) * t40;
t35 = pkin(2) * t49 + pkin(1);
t39 = sin(pkin(10));
t41 = cos(pkin(10));
t60 = t41 * t25 + t39 * t35 + r_base(2);
t59 = t42 * pkin(7) + t63;
t26 = -t49 * t38 - t46 * t64;
t24 = t26 * t42;
t14 = -t24 * t41 + t27 * t39;
t58 = t14 * pkin(3) + t60;
t57 = -t25 * t39 + t41 * t35 + r_base(1);
t56 = t40 * t46 * pkin(2) + t42 * qJ(3) + t59;
t16 = t24 * t39 + t27 * t41;
t55 = t16 * pkin(3) + t57;
t23 = t26 * t40;
t54 = -t23 * pkin(3) + t56;
t52 = t42 * t27;
t13 = t39 * t26 + t41 * t52;
t53 = -pkin(8) * t13 + t58;
t15 = t26 * t41 - t39 * t52;
t51 = -pkin(8) * t15 + t55;
t22 = t27 * t40;
t50 = -t22 * pkin(8) + t54;
t47 = cos(qJ(5));
t34 = pkin(5) * t47 + pkin(4);
t18 = -t23 * t48 + t42 * t45;
t10 = t16 * t48 + t39 * t69;
t8 = t14 * t48 - t41 * t69;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t63 - mrSges(2,3) - m(3) * t59 - (t46 * mrSges(3,1) + t83) * t40 - m(4) * t56 + t23 * mrSges(4,1) - m(5) * t50 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t50) - m(7) * (t18 * t34 + t54) + t80 * (t18 * t47 - t22 * t44) + t79 * (-t18 * t44 - t22 * t47) + t78 * t42 + t75 * t22 + t73 * (-t23 * t45 - t42 * t48)) * g(3) + (-m(7) * (t34 * t8 + t58) - m(6) * (pkin(4) * t8 + t53) - m(4) * t60 - m(5) * t53 - mrSges(1,2) - t8 * mrSges(5,1) + t72 * t39 + t77 * r_base(2) + t80 * (-t13 * t44 + t47 * t8) + t75 * t13 + t79 * (-t13 * t47 - t44 * t8) + t82 * t14 + t81 * t41) * g(2) + (-m(7) * (t10 * t34 + t55) - m(6) * (pkin(4) * t10 + t51) - m(5) * t51 - m(4) * t57 - mrSges(1,1) - t10 * mrSges(5,1) + t72 * t41 + t77 * r_base(1) + t80 * (t10 * t47 - t15 * t44) + t79 * (-t10 * t44 - t15 * t47) + t75 * t15 + t82 * t16 - t81 * t39) * g(1);
U  = t1;
