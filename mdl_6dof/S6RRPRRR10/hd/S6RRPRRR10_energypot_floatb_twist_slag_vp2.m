% Calculate potential energy for
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:03
% EndTime: 2019-03-09 14:18:04
% DurationCPUTime: 0.94s
% Computational Cost: add. (376->115), mult. (559->126), div. (0->0), fcn. (633->14), ass. (0->46)
t33 = qJ(5) + qJ(6);
t27 = sin(t33);
t28 = cos(t33);
t38 = -pkin(9) - qJ(3);
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t71 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t72 = -m(5) - m(6) - m(7);
t74 = -t27 * mrSges(7,1) - t42 * mrSges(6,2) - t28 * mrSges(7,2) - t72 * t38 - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t39 - t71;
t73 = -m(1) - m(2);
t69 = -m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t24 = pkin(5) * t42 + pkin(4);
t67 = -m(6) * pkin(4) - m(7) * t24 - t42 * mrSges(6,1) - mrSges(7,1) * t28 + t39 * mrSges(6,2) + mrSges(7,2) * t27 - mrSges(5,1);
t34 = sin(pkin(12));
t66 = pkin(3) * t34;
t35 = sin(pkin(6));
t40 = sin(qJ(2));
t63 = t35 * t40;
t41 = sin(qJ(1));
t62 = t35 * t41;
t43 = cos(qJ(2));
t61 = t35 * t43;
t44 = cos(qJ(1));
t60 = t35 * t44;
t59 = t40 * t41;
t58 = t40 * t44;
t57 = t41 * t43;
t56 = t43 * t44;
t55 = pkin(7) + r_base(3);
t54 = t41 * pkin(1) + r_base(2);
t53 = t34 * t62;
t52 = t39 * t61;
t37 = cos(pkin(6));
t51 = t37 * pkin(8) + t55;
t50 = t44 * pkin(1) + pkin(8) * t62 + r_base(1);
t49 = -pkin(8) * t60 + t54;
t36 = cos(pkin(12));
t23 = pkin(3) * t36 + pkin(2);
t48 = t23 * t63 + t37 * t66 + t38 * t61 + t51;
t32 = pkin(12) + qJ(4);
t26 = cos(t32);
t25 = sin(t32);
t14 = -t37 * t59 + t56;
t12 = t37 * t58 + t57;
t8 = t25 * t37 + t26 * t63;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t55 - mrSges(2,3) - m(5) * t48 - t8 * mrSges(5,1) + mrSges(5,3) * t61 - m(6) * (pkin(4) * t8 + t48) - (t42 * t8 - t52) * mrSges(6,1) - (-t39 * t8 - t42 * t61) * mrSges(6,2) - m(7) * (-pkin(5) * t52 + t8 * t24 + t48) - (-t27 * t61 + t28 * t8) * mrSges(7,1) - (-t27 * t8 - t28 * t61) * mrSges(7,2) + (-m(3) - m(4)) * t51 + (-t34 * mrSges(4,1) - t36 * mrSges(4,2) - mrSges(3,3)) * t37 + (t71 * t43 + (-m(4) * pkin(2) - t36 * mrSges(4,1) + t34 * mrSges(4,2) - mrSges(3,1)) * t40) * t35 + t69 * (t25 * t63 - t37 * t26)) * g(3) + (-mrSges(1,2) - m(4) * (pkin(2) * t12 + t49) - m(3) * t49 + mrSges(3,3) * t60 - (-t12 * t34 - t36 * t60) * mrSges(4,2) - (t12 * t36 - t34 * t60) * mrSges(4,1) - t12 * mrSges(3,1) - t41 * mrSges(2,1) - t44 * mrSges(2,2) + t73 * r_base(2) + t72 * (t12 * t23 + (-pkin(8) - t66) * t60 + t54) + t67 * (t12 * t26 - t25 * t60) + t69 * (t12 * t25 + t26 * t60) + t74 * (-t37 * t56 + t59)) * g(2) + (-mrSges(3,3) * t62 - mrSges(1,1) - m(4) * (pkin(2) * t14 + t50) - m(3) * t50 - (t14 * t36 + t53) * mrSges(4,1) - (-t14 * t34 + t36 * t62) * mrSges(4,2) - t14 * mrSges(3,1) + t41 * mrSges(2,2) - t44 * mrSges(2,1) + t73 * r_base(1) + t72 * (pkin(3) * t53 + t14 * t23 + t50) + t69 * (t14 * t25 - t26 * t62) + t67 * (t14 * t26 + t25 * t62) + t74 * (t37 * t57 + t58)) * g(1);
U  = t1;
