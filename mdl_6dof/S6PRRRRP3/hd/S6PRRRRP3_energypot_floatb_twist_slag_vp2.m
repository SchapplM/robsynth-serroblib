% Calculate potential energy for
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:15
% EndTime: 2019-03-09 00:06:15
% DurationCPUTime: 0.72s
% Computational Cost: add. (344->109), mult. (618->123), div. (0->0), fcn. (717->12), ass. (0->50)
t74 = -m(2) - m(1);
t73 = -m(4) - m(5);
t72 = -m(6) - m(7);
t38 = sin(qJ(4));
t56 = pkin(4) * t38 + pkin(8);
t71 = -mrSges(6,1) - mrSges(7,1);
t70 = -mrSges(6,2) - mrSges(7,2);
t43 = -pkin(10) - pkin(9);
t69 = -m(5) * pkin(9) + m(6) * t43 + m(7) * (-qJ(6) + t43) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t41 = cos(qJ(4));
t34 = qJ(4) + qJ(5);
t27 = sin(t34);
t65 = pkin(5) * t27 + t56;
t68 = -m(6) * t56 - m(7) * t65 - t38 * mrSges(5,1) - t41 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t26 = t41 * pkin(4) + pkin(3);
t28 = cos(t34);
t19 = pkin(5) * t28 + t26;
t67 = -m(5) * pkin(3) - m(6) * t26 - m(7) * t19 - t41 * mrSges(5,1) + t38 * mrSges(5,2) - mrSges(4,1);
t64 = cos(qJ(3));
t35 = sin(pkin(11));
t36 = sin(pkin(6));
t63 = t35 * t36;
t37 = cos(pkin(11));
t62 = t36 * t37;
t39 = sin(qJ(3));
t61 = t36 * t39;
t40 = sin(qJ(2));
t60 = t36 * t40;
t42 = cos(qJ(2));
t59 = t36 * t42;
t58 = cos(pkin(6));
t57 = qJ(1) + r_base(3);
t55 = t36 * t64;
t54 = t40 * t58;
t53 = t42 * t58;
t52 = t37 * pkin(1) + pkin(7) * t63 + r_base(1);
t51 = t58 * pkin(7) + t57;
t16 = -t35 * t54 + t37 * t42;
t50 = t16 * pkin(2) + t52;
t49 = pkin(2) * t60 + t51;
t48 = t35 * pkin(1) - pkin(7) * t62 + r_base(2);
t14 = t35 * t42 + t37 * t54;
t47 = t14 * pkin(2) + t48;
t45 = -pkin(8) * t59 + t49;
t18 = t58 * t39 + t40 * t55;
t15 = t35 * t53 + t37 * t40;
t13 = t35 * t40 - t37 * t53;
t10 = t16 * t64 + t35 * t61;
t8 = t14 * t64 - t37 * t61;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t57 - mrSges(2,3) - m(3) * t51 - t58 * mrSges(3,3) - (t40 * mrSges(3,1) + t42 * mrSges(3,2)) * t36 - m(4) * t45 - t18 * mrSges(4,1) + mrSges(4,3) * t59 - m(5) * (pkin(3) * t18 + t45) - (t18 * t41 - t38 * t59) * mrSges(5,1) - (-t18 * t38 - t41 * t59) * mrSges(5,2) - m(6) * (t18 * t26 - t56 * t59 + t49) - m(7) * (t18 * t19 - t65 * t59 + t49) + t71 * (t18 * t28 - t27 * t59) + t70 * (-t18 * t27 - t28 * t59) + t69 * (t39 * t60 - t58 * t64)) * g(3) + (-m(3) * t48 - t35 * mrSges(2,1) - t14 * mrSges(3,1) - t37 * mrSges(2,2) + mrSges(3,3) * t62 - mrSges(1,2) + t74 * r_base(2) + t72 * t47 + t73 * (t13 * pkin(8) + t47) + t71 * (t13 * t27 + t28 * t8) + t67 * t8 + t68 * t13 + t70 * (t13 * t28 - t27 * t8) + t69 * (t14 * t39 + t37 * t55)) * g(2) + (-m(3) * t52 - t37 * mrSges(2,1) - t16 * mrSges(3,1) + t35 * mrSges(2,2) - mrSges(3,3) * t63 - mrSges(1,1) + t74 * r_base(1) + t72 * t50 + t73 * (t15 * pkin(8) + t50) + t71 * (t10 * t28 + t15 * t27) + t70 * (-t10 * t27 + t15 * t28) + t67 * t10 + t68 * t15 + t69 * (t16 * t39 - t35 * t55)) * g(1);
U  = t1;
