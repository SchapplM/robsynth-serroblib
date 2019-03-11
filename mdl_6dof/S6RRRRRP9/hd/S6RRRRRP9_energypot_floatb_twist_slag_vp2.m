% Calculate potential energy for
% S6RRRRRP9
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:21
% EndTime: 2019-03-10 02:00:21
% DurationCPUTime: 0.68s
% Computational Cost: add. (344->109), mult. (618->122), div. (0->0), fcn. (717->12), ass. (0->49)
t73 = -m(1) - m(2);
t72 = -m(4) - m(5);
t71 = -m(6) - m(7);
t36 = sin(qJ(4));
t56 = pkin(4) * t36 + pkin(9);
t70 = -mrSges(6,1) - mrSges(7,1);
t69 = -mrSges(6,2) - mrSges(7,2);
t43 = -pkin(11) - pkin(10);
t68 = -m(5) * pkin(10) + m(6) * t43 + m(7) * (-qJ(6) + t43) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t40 = cos(qJ(4));
t34 = qJ(4) + qJ(5);
t27 = sin(t34);
t64 = pkin(5) * t27 + t56;
t67 = -m(6) * t56 - m(7) * t64 - t36 * mrSges(5,1) - t40 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t26 = t40 * pkin(4) + pkin(3);
t28 = cos(t34);
t19 = pkin(5) * t28 + t26;
t66 = -m(5) * pkin(3) - m(6) * t26 - m(7) * t19 - t40 * mrSges(5,1) + t36 * mrSges(5,2) - mrSges(4,1);
t63 = cos(qJ(3));
t35 = sin(pkin(6));
t38 = sin(qJ(2));
t62 = t35 * t38;
t39 = sin(qJ(1));
t61 = t35 * t39;
t41 = cos(qJ(2));
t60 = t35 * t41;
t42 = cos(qJ(1));
t59 = t35 * t42;
t58 = cos(pkin(6));
t57 = pkin(7) + r_base(3);
t55 = t35 * t63;
t54 = t58 * pkin(8) + t57;
t53 = t39 * t58;
t52 = t42 * t58;
t51 = t42 * pkin(1) + pkin(8) * t61 + r_base(1);
t50 = pkin(2) * t62 + t54;
t18 = -t38 * t53 + t42 * t41;
t49 = t18 * pkin(2) + t51;
t48 = t39 * pkin(1) - pkin(8) * t59 + r_base(2);
t16 = t38 * t52 + t39 * t41;
t47 = t16 * pkin(2) + t48;
t45 = -pkin(9) * t60 + t50;
t37 = sin(qJ(3));
t17 = t42 * t38 + t41 * t53;
t15 = t38 * t39 - t41 * t52;
t14 = t37 * t58 + t38 * t55;
t10 = t18 * t63 + t37 * t61;
t8 = t16 * t63 - t37 * t59;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t57 - mrSges(2,3) - m(3) * t54 - t58 * mrSges(3,3) - (t38 * mrSges(3,1) + t41 * mrSges(3,2)) * t35 - m(4) * t45 - t14 * mrSges(4,1) + mrSges(4,3) * t60 - m(5) * (pkin(3) * t14 + t45) - (t14 * t40 - t36 * t60) * mrSges(5,1) - (-t14 * t36 - t40 * t60) * mrSges(5,2) - m(6) * (t14 * t26 - t56 * t60 + t50) - m(7) * (t14 * t19 - t64 * t60 + t50) + t70 * (t14 * t28 - t27 * t60) + t69 * (-t14 * t27 - t28 * t60) + t68 * (t37 * t62 - t58 * t63)) * g(3) + (-m(3) * t48 - t39 * mrSges(2,1) - t16 * mrSges(3,1) - t42 * mrSges(2,2) + mrSges(3,3) * t59 - mrSges(1,2) + t73 * r_base(2) + t71 * t47 + t72 * (t15 * pkin(9) + t47) + t70 * (t15 * t27 + t28 * t8) + t66 * t8 + t67 * t15 + t69 * (t15 * t28 - t27 * t8) + t68 * (t16 * t37 + t42 * t55)) * g(2) + (-m(3) * t51 - t42 * mrSges(2,1) - t18 * mrSges(3,1) + t39 * mrSges(2,2) - mrSges(3,3) * t61 - mrSges(1,1) + t73 * r_base(1) + t71 * t49 + t72 * (t17 * pkin(9) + t49) + t70 * (t10 * t28 + t17 * t27) + t69 * (-t10 * t27 + t17 * t28) + t66 * t10 + t67 * t17 + t68 * (t18 * t37 - t39 * t55)) * g(1);
U  = t1;
