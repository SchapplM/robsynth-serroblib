% Calculate potential energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:28
% EndTime: 2019-03-09 12:55:28
% DurationCPUTime: 0.70s
% Computational Cost: add. (281->108), mult. (550->118), div. (0->0), fcn. (622->10), ass. (0->51)
t74 = -m(1) - m(2);
t73 = -mrSges(3,1) + mrSges(4,2);
t72 = -mrSges(6,1) - mrSges(7,1);
t71 = -mrSges(6,2) - mrSges(7,2);
t70 = mrSges(3,3) + mrSges(4,1);
t69 = -mrSges(4,3) + mrSges(3,2);
t34 = sin(qJ(5));
t68 = -m(7) * (pkin(5) * t34 + pkin(9)) - mrSges(5,3) + t73;
t67 = m(6) * pkin(10) - m(7) * (-qJ(6) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t31 = sin(pkin(6));
t36 = sin(qJ(2));
t66 = t31 * t36;
t37 = sin(qJ(1));
t65 = t31 * t37;
t40 = cos(qJ(2));
t64 = t31 * t40;
t41 = cos(qJ(1));
t63 = t31 * t41;
t62 = t34 * t36;
t61 = t36 * t37;
t60 = t36 * t41;
t59 = t37 * t40;
t58 = t40 * t41;
t57 = qJ(3) * t40;
t56 = pkin(7) + r_base(3);
t55 = pkin(8) * t63;
t54 = t37 * pkin(1) + r_base(2);
t32 = cos(pkin(6));
t52 = t32 * pkin(8) + t56;
t51 = t41 * pkin(1) + pkin(8) * t65 + r_base(1);
t50 = pkin(2) * t66 + t52;
t49 = t32 * pkin(3) + pkin(9) * t66 + t50;
t16 = -t32 * t58 + t61;
t17 = t32 * t60 + t59;
t48 = t17 * pkin(2) + t16 * qJ(3) + t54;
t18 = t32 * t59 + t60;
t19 = -t32 * t61 + t58;
t47 = t19 * pkin(2) + qJ(3) * t18 + t51;
t46 = pkin(3) * t65 + t47;
t45 = -t31 * t57 + t49;
t44 = pkin(9) * t19 + t46;
t43 = (-pkin(3) - pkin(8)) * t63 + t48;
t42 = t17 * pkin(9) + t43;
t39 = cos(qJ(4));
t38 = cos(qJ(5));
t35 = sin(qJ(4));
t26 = pkin(5) * t38 + pkin(4);
t15 = t32 * t39 - t35 * t64;
t10 = t16 * t35 - t39 * t63;
t8 = t18 * t35 + t39 * t65;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t56 - mrSges(2,3) - m(3) * t52 - m(4) * t50 - m(5) * t45 - t15 * mrSges(5,1) - mrSges(5,3) * t66 - m(6) * (pkin(4) * t15 + t45) - m(7) * (t15 * t26 + t49) + t72 * (t15 * t38 + t31 * t62) + t71 * (-t15 * t34 + t38 * t66) - t70 * t32 + (-m(7) * (pkin(5) * t62 - t57) + (m(4) * qJ(3) - t69) * t40 + t73 * t36) * t31 - t67 * (t32 * t35 + t39 * t64)) * g(3) + (-m(3) * (t54 - t55) - m(4) * (t48 - t55) - m(7) * (t10 * t26 + t43) - m(6) * (t10 * pkin(4) + t42) - m(5) * t42 - mrSges(1,2) - t10 * mrSges(5,1) - t37 * mrSges(2,1) - t41 * mrSges(2,2) + t74 * r_base(2) + t70 * t63 + t72 * (t10 * t38 + t17 * t34) + t71 * (-t10 * t34 + t17 * t38) + t69 * t16 + t67 * (t16 * t39 + t35 * t63) + t68 * t17) * g(2) + (-m(7) * (t26 * t8 + t46) - m(6) * (pkin(4) * t8 + t44) - m(3) * t51 - m(4) * t47 - m(5) * t44 - mrSges(1,1) - t8 * mrSges(5,1) + t37 * mrSges(2,2) - t41 * mrSges(2,1) + t74 * r_base(1) - t70 * t65 + t72 * (t19 * t34 + t38 * t8) + t69 * t18 + t71 * (t19 * t38 - t34 * t8) - t67 * (-t18 * t39 + t35 * t65) + t68 * t19) * g(1);
U  = t1;
