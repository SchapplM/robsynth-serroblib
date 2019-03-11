% Calculate potential energy for
% S6PRPRRP3
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:35
% EndTime: 2019-03-08 20:04:36
% DurationCPUTime: 0.81s
% Computational Cost: add. (364->119), mult. (559->134), div. (0->0), fcn. (633->12), ass. (0->49)
t73 = -m(1) - m(2);
t72 = -mrSges(6,1) - mrSges(7,1);
t71 = -mrSges(6,2) - mrSges(7,2);
t70 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t69 = -mrSges(5,3) - t70;
t68 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t37 = sin(pkin(11));
t67 = pkin(3) * t37;
t38 = sin(pkin(10));
t41 = cos(pkin(10));
t46 = sin(qJ(2));
t42 = cos(pkin(6));
t48 = cos(qJ(2));
t59 = t42 * t48;
t17 = t38 * t46 - t41 * t59;
t45 = sin(qJ(5));
t66 = t17 * t45;
t19 = t38 * t59 + t41 * t46;
t65 = t19 * t45;
t39 = sin(pkin(6));
t64 = t38 * t39;
t63 = t39 * t41;
t62 = t39 * t46;
t61 = t39 * t48;
t60 = t42 * t46;
t58 = pkin(1) * t38 + r_base(2);
t57 = t37 * t64;
t56 = t45 * t61;
t55 = qJ(1) + r_base(3);
t54 = pkin(1) * t41 + pkin(7) * t64 + r_base(1);
t53 = pkin(7) * t42 + t55;
t52 = -pkin(7) * t63 + t58;
t20 = -t38 * t60 + t41 * t48;
t40 = cos(pkin(11));
t29 = pkin(3) * t40 + pkin(2);
t44 = -pkin(8) - qJ(3);
t51 = pkin(3) * t57 - t19 * t44 + t20 * t29 + t54;
t50 = t29 * t62 + t42 * t67 + t44 * t61 + t53;
t18 = t38 * t48 + t41 * t60;
t49 = t18 * t29 - t17 * t44 + (-pkin(7) - t67) * t63 + t58;
t47 = cos(qJ(5));
t36 = pkin(11) + qJ(4);
t32 = cos(t36);
t31 = sin(t36);
t30 = pkin(5) * t47 + pkin(4);
t14 = t31 * t42 + t32 * t62;
t8 = t20 * t32 + t31 * t64;
t6 = t18 * t32 - t31 * t63;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t55 - mrSges(2,3) - m(5) * t50 - t14 * mrSges(5,1) + mrSges(5,3) * t61 - m(6) * (pkin(4) * t14 + t50) - m(7) * (-pkin(5) * t56 + t14 * t30 + t50) + t71 * (-t14 * t45 - t47 * t61) + (-m(3) - m(4)) * t53 + (-t37 * mrSges(4,1) - t40 * mrSges(4,2) - mrSges(3,3)) * t42 + (t70 * t48 + (-m(4) * pkin(2) - mrSges(4,1) * t40 + mrSges(4,2) * t37 - mrSges(3,1)) * t46) * t39 + t72 * (t14 * t47 - t56) + t68 * (t31 * t62 - t32 * t42)) * g(3) + (-mrSges(1,2) - t6 * mrSges(5,1) - t18 * mrSges(3,1) - t38 * mrSges(2,1) - t41 * mrSges(2,2) - m(6) * (pkin(4) * t6 + t49) - m(5) * t49 - m(4) * (pkin(2) * t18 + t52) - m(3) * t52 + mrSges(3,3) * t63 - (-t18 * t37 - t40 * t63) * mrSges(4,2) - (t18 * t40 - t37 * t63) * mrSges(4,1) - m(7) * (pkin(5) * t66 + t30 * t6 + t49) + t73 * r_base(2) + t72 * (t47 * t6 + t66) + t71 * (t17 * t47 - t45 * t6) + t68 * (t18 * t31 + t32 * t63) + t69 * t17) * g(2) + (-mrSges(3,3) * t64 - mrSges(1,1) - t8 * mrSges(5,1) - t20 * mrSges(3,1) + t38 * mrSges(2,2) - t41 * mrSges(2,1) - m(6) * (pkin(4) * t8 + t51) - m(5) * t51 - m(4) * (pkin(2) * t20 + t54) - m(3) * t54 - (t20 * t40 + t57) * mrSges(4,1) - (-t20 * t37 + t40 * t64) * mrSges(4,2) - m(7) * (pkin(5) * t65 + t30 * t8 + t51) + t73 * r_base(1) + t72 * (t47 * t8 + t65) + t71 * (t19 * t47 - t45 * t8) + t68 * (t20 * t31 - t32 * t64) + t69 * t19) * g(1);
U  = t1;
