% Calculate potential energy for
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:02
% EndTime: 2019-03-08 20:18:02
% DurationCPUTime: 0.66s
% Computational Cost: add. (296->99), mult. (596->112), div. (0->0), fcn. (686->10), ass. (0->49)
t76 = -m(1) - m(2);
t75 = -m(6) - m(7);
t74 = -mrSges(3,1) + mrSges(4,2);
t73 = mrSges(3,3) + mrSges(4,1);
t72 = -mrSges(4,3) + mrSges(3,2);
t71 = -mrSges(5,3) + t74;
t70 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t69 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t68 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t36 = sin(pkin(10));
t37 = sin(pkin(6));
t67 = t36 * t37;
t38 = cos(pkin(10));
t66 = t37 * t38;
t41 = sin(qJ(4));
t65 = t37 * t41;
t42 = sin(qJ(2));
t64 = t37 * t42;
t44 = cos(qJ(4));
t63 = t37 * t44;
t45 = cos(qJ(2));
t62 = t37 * t45;
t39 = cos(pkin(6));
t61 = t39 * t42;
t60 = t39 * t45;
t59 = pkin(7) * t66;
t58 = t36 * pkin(1) + r_base(2);
t57 = qJ(1) + r_base(3);
t56 = t38 * pkin(1) + pkin(7) * t67 + r_base(1);
t55 = t39 * pkin(7) + t57;
t54 = pkin(2) * t64 + t55;
t19 = t36 * t42 - t38 * t60;
t20 = t36 * t45 + t38 * t61;
t53 = t20 * pkin(2) + qJ(3) * t19 + t58;
t21 = t36 * t60 + t38 * t42;
t22 = -t36 * t61 + t38 * t45;
t52 = t22 * pkin(2) + qJ(3) * t21 + t56;
t51 = t39 * pkin(3) + pkin(8) * t64 - qJ(3) * t62 + t54;
t50 = pkin(3) * t67 + pkin(8) * t22 + t52;
t49 = pkin(8) * t20 + (-pkin(3) - pkin(7)) * t66 + t53;
t43 = cos(qJ(5));
t40 = sin(qJ(5));
t24 = t39 * t44 - t41 * t62;
t23 = t39 * t41 + t44 * t62;
t10 = t19 * t41 - t38 * t63;
t9 = t19 * t44 + t38 * t65;
t8 = t21 * t41 + t36 * t63;
t7 = -t21 * t44 + t36 * t65;
t1 = (-m(1) * r_base(3) - m(2) * t57 - m(3) * t55 - m(4) * t54 - m(5) * t51 - t24 * mrSges(5,1) - mrSges(5,3) * t64 - mrSges(1,3) - mrSges(2,3) + t75 * (t24 * pkin(4) + t23 * pkin(9) + t51) - t73 * t39 + ((m(4) * qJ(3) - t72) * t45 + t74 * t42) * t37 + t69 * (t24 * t43 + t40 * t64) + t68 * (t24 * t40 - t43 * t64) + t70 * t23) * g(3) + (-m(3) * (t58 - t59) - m(4) * (t53 - t59) - m(5) * t49 - mrSges(1,2) - t10 * mrSges(5,1) - t36 * mrSges(2,1) - t38 * mrSges(2,2) + t76 * r_base(2) + t73 * t66 + t75 * (t10 * pkin(4) - pkin(9) * t9 + t49) + t69 * (t10 * t43 + t20 * t40) + t68 * (t10 * t40 - t20 * t43) + t72 * t19 - t70 * t9 + t71 * t20) * g(2) + (-m(3) * t56 - m(4) * t52 - m(5) * t50 - t38 * mrSges(2,1) - t8 * mrSges(5,1) + t36 * mrSges(2,2) - mrSges(1,1) + t76 * r_base(1) - t73 * t67 + t75 * (t8 * pkin(4) + pkin(9) * t7 + t50) + t72 * t21 + t69 * (t22 * t40 + t43 * t8) + t68 * (-t22 * t43 + t40 * t8) + t70 * t7 + t71 * t22) * g(1);
U  = t1;
