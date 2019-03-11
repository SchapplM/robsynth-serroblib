% Calculate potential energy for
% S6RRPRRR5
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:37:55
% EndTime: 2019-03-09 13:37:56
% DurationCPUTime: 0.80s
% Computational Cost: add. (421->100), mult. (867->115), div. (0->0), fcn. (1059->14), ass. (0->51)
t77 = -m(5) - m(6);
t76 = -mrSges(3,3) - mrSges(4,3);
t35 = sin(pkin(12));
t40 = sin(qJ(2));
t44 = cos(qJ(2));
t61 = cos(pkin(12));
t21 = -t40 * t35 + t44 * t61;
t75 = -m(3) - m(1) - m(2);
t74 = -m(3) * pkin(1) - mrSges(2,1);
t73 = m(3) * pkin(8) - t76;
t72 = -m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t34 = qJ(5) + qJ(6);
t31 = sin(t34);
t32 = cos(t34);
t38 = sin(qJ(5));
t42 = cos(qJ(5));
t71 = m(7) * (pkin(5) * t38 + pkin(9)) + t38 * mrSges(6,1) + t31 * mrSges(7,1) + t42 * mrSges(6,2) + t32 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t70 = -m(6) * pkin(4) - m(7) * (pkin(5) * t42 + pkin(4)) - t42 * mrSges(6,1) - t32 * mrSges(7,1) + t38 * mrSges(6,2) + t31 * mrSges(7,2) - mrSges(5,1);
t69 = pkin(2) * t40;
t36 = sin(pkin(6));
t41 = sin(qJ(1));
t68 = t36 * t41;
t45 = cos(qJ(1));
t67 = t36 * t45;
t65 = t40 * t41;
t64 = t40 * t45;
t63 = t41 * t44;
t62 = t44 * t45;
t60 = pkin(7) + r_base(3);
t37 = cos(pkin(6));
t58 = t37 * pkin(8) + t60;
t19 = t37 * t69 + (-pkin(8) - qJ(3)) * t36;
t29 = pkin(2) * t44 + pkin(1);
t56 = t45 * t19 + t41 * t29 + r_base(2);
t20 = -t44 * t35 - t40 * t61;
t18 = t20 * t37;
t8 = -t18 * t45 + t21 * t41;
t55 = t8 * pkin(3) + t56;
t54 = -t41 * t19 + t45 * t29 + r_base(1);
t53 = t37 * qJ(3) + t36 * t69 + t58;
t10 = t18 * t41 + t21 * t45;
t52 = t10 * pkin(3) + t54;
t17 = t20 * t36;
t51 = -t17 * pkin(3) + t53;
t48 = t37 * t21;
t43 = cos(qJ(4));
t39 = sin(qJ(4));
t16 = t21 * t36;
t9 = t20 * t45 - t41 * t48;
t7 = t41 * t20 + t45 * t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t60 - mrSges(2,3) - m(3) * t58 - (t40 * mrSges(3,1) + t44 * mrSges(3,2)) * t36 - m(4) * t53 + t17 * mrSges(4,1) - m(7) * t51 + t77 * (-t16 * pkin(9) + t51) + t76 * t37 + t70 * (-t17 * t43 + t37 * t39) + t71 * t16 + t72 * (-t17 * t39 - t37 * t43)) * g(3) + (-(t37 * t64 + t63) * mrSges(3,1) - (t37 * t62 - t65) * mrSges(3,2) - m(7) * t55 - m(4) * t56 - mrSges(1,2) - t8 * mrSges(4,1) - t45 * mrSges(2,2) + t74 * t41 + t75 * r_base(2) + t77 * (-t7 * pkin(9) + t55) + t70 * (-t39 * t67 + t43 * t8) + t71 * t7 + t73 * t67 + t72 * (t39 * t8 + t43 * t67)) * g(2) + (-(-t37 * t63 - t64) * mrSges(3,2) - (-t37 * t65 + t62) * mrSges(3,1) - m(7) * t52 - m(4) * t54 - mrSges(1,1) - t10 * mrSges(4,1) + t41 * mrSges(2,2) + t74 * t45 + t75 * r_base(1) + t77 * (-t9 * pkin(9) + t52) + t70 * (t10 * t43 + t39 * t68) + t71 * t9 - t73 * t68 + t72 * (t10 * t39 - t43 * t68)) * g(1);
U  = t1;
