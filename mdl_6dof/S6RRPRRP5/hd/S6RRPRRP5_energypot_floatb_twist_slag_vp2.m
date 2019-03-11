% Calculate potential energy for
% S6RRPRRP5
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:40
% EndTime: 2019-03-09 11:57:41
% DurationCPUTime: 0.79s
% Computational Cost: add. (409->111), mult. (867->132), div. (0->0), fcn. (1059->12), ass. (0->55)
t80 = -mrSges(6,1) - mrSges(7,1);
t79 = -mrSges(6,2) - mrSges(7,2);
t78 = -mrSges(3,3) - mrSges(4,3);
t38 = sin(pkin(11));
t44 = sin(qJ(2));
t48 = cos(qJ(2));
t64 = cos(pkin(11));
t27 = -t44 * t38 + t48 * t64;
t77 = -m(3) - m(2) - m(1);
t76 = -m(3) * pkin(1) - mrSges(2,1);
t75 = m(3) * pkin(8) - t78;
t42 = sin(qJ(5));
t74 = m(7) * (pkin(5) * t42 + pkin(9)) - mrSges(4,2) + mrSges(5,3);
t73 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t72 = pkin(2) * t44;
t39 = sin(pkin(6));
t45 = sin(qJ(1));
t71 = t39 * t45;
t49 = cos(qJ(1));
t70 = t39 * t49;
t68 = t44 * t49;
t67 = t45 * t44;
t66 = t45 * t48;
t65 = t48 * t49;
t63 = pkin(7) + r_base(3);
t40 = cos(pkin(6));
t61 = t40 * pkin(8) + t63;
t25 = t40 * t72 + (-pkin(8) - qJ(3)) * t39;
t35 = pkin(2) * t48 + pkin(1);
t59 = t49 * t25 + t45 * t35 + r_base(2);
t26 = -t48 * t38 - t44 * t64;
t24 = t26 * t40;
t14 = -t24 * t49 + t45 * t27;
t58 = t14 * pkin(3) + t59;
t57 = -t25 * t45 + t49 * t35 + r_base(1);
t56 = t40 * qJ(3) + t39 * t72 + t61;
t16 = t45 * t24 + t27 * t49;
t55 = t16 * pkin(3) + t57;
t23 = t26 * t39;
t54 = -t23 * pkin(3) + t56;
t52 = t40 * t27;
t13 = t45 * t26 + t49 * t52;
t53 = -pkin(9) * t13 + t58;
t15 = t26 * t49 - t45 * t52;
t51 = -pkin(9) * t15 + t55;
t22 = t27 * t39;
t50 = -pkin(9) * t22 + t54;
t47 = cos(qJ(4));
t46 = cos(qJ(5));
t43 = sin(qJ(4));
t34 = pkin(5) * t46 + pkin(4);
t18 = -t23 * t47 + t40 * t43;
t10 = t16 * t47 + t43 * t71;
t8 = t14 * t47 - t43 * t70;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t63 - mrSges(2,3) - m(3) * t61 - (t44 * mrSges(3,1) + t48 * mrSges(3,2)) * t39 - m(4) * t56 + t23 * mrSges(4,1) - m(5) * t50 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t50) - m(7) * (t18 * t34 + t54) + t80 * (t18 * t46 - t22 * t42) + t79 * (-t18 * t42 - t22 * t46) + t78 * t40 + t74 * t22 + t73 * (-t23 * t43 - t40 * t47)) * g(3) + (-(t40 * t65 - t67) * mrSges(3,2) - (t40 * t68 + t66) * mrSges(3,1) - m(7) * (t34 * t8 + t58) - m(4) * t59 - m(6) * (pkin(4) * t8 + t53) - m(5) * t53 - t49 * mrSges(2,2) - mrSges(1,2) - t8 * mrSges(5,1) - t14 * mrSges(4,1) + t76 * t45 + t77 * r_base(2) + t75 * t70 + t80 * (-t13 * t42 + t46 * t8) + t74 * t13 + t79 * (-t13 * t46 - t42 * t8) + t73 * (t14 * t43 + t47 * t70)) * g(2) + (-(-t40 * t67 + t65) * mrSges(3,1) - (-t40 * t66 - t68) * mrSges(3,2) - m(7) * (t10 * t34 + t55) - m(4) * t57 - m(6) * (pkin(4) * t10 + t51) - m(5) * t51 + t45 * mrSges(2,2) - mrSges(1,1) - t10 * mrSges(5,1) - t16 * mrSges(4,1) + t76 * t49 + t77 * r_base(1) - t75 * t71 + t80 * (t10 * t46 - t15 * t42) + t79 * (-t10 * t42 - t15 * t46) + t74 * t15 + t73 * (t16 * t43 - t47 * t71)) * g(1);
U  = t1;
