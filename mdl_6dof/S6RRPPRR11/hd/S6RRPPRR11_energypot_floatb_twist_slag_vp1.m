% Calculate potential energy for
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:38:55
% EndTime: 2019-03-09 09:38:56
% DurationCPUTime: 0.73s
% Computational Cost: add. (305->129), mult. (482->150), div. (0->0), fcn. (547->12), ass. (0->54)
t69 = pkin(10) + rSges(7,3);
t68 = qJ(4) + rSges(5,3);
t31 = sin(pkin(11));
t67 = pkin(4) * t31;
t38 = sin(qJ(1));
t66 = g(1) * t38;
t41 = cos(qJ(1));
t65 = g(2) * t41;
t32 = sin(pkin(6));
t64 = t32 * t38;
t40 = cos(qJ(2));
t63 = t32 * t40;
t62 = t32 * t41;
t37 = sin(qJ(2));
t61 = t37 * t38;
t60 = t37 * t41;
t59 = t38 * t40;
t58 = t40 * t41;
t57 = pkin(7) + r_base(3);
t56 = t38 * pkin(1) + r_base(2);
t34 = cos(pkin(6));
t55 = t34 * pkin(8) + t57;
t15 = t34 * t60 + t59;
t54 = t15 * pkin(2) + t56;
t53 = t41 * pkin(1) + pkin(8) * t64 + r_base(1);
t33 = cos(pkin(11));
t24 = pkin(4) * t33 + pkin(3);
t52 = (-pkin(8) - t24) * t65;
t51 = t32 * t37 * pkin(2) + t55;
t17 = -t34 * t61 + t58;
t50 = t17 * pkin(2) + t53;
t49 = (-qJ(3) - t67) * t40;
t48 = t34 * t24 + t51;
t47 = rSges(5,1) * t33 - rSges(5,2) * t31 + pkin(3);
t46 = rSges(5,1) * t31 + rSges(5,2) * t33 + qJ(3);
t14 = -t34 * t58 + t61;
t45 = t14 * qJ(3) + t54;
t16 = t34 * t59 + t60;
t44 = t16 * qJ(3) + t50;
t35 = -pkin(9) - qJ(4);
t43 = t14 * t67 - t15 * t35 + t45;
t42 = t16 * t67 - t17 * t35 + t24 * t64 + t44;
t39 = cos(qJ(6));
t36 = sin(qJ(6));
t30 = pkin(11) + qJ(5);
t26 = cos(t30);
t25 = sin(t30);
t7 = -t25 * t63 + t26 * t34;
t6 = t25 * t34 + t26 * t63;
t4 = t14 * t25 - t26 * t62;
t3 = t14 * t26 + t25 * t62;
t2 = t16 * t25 + t26 * t64;
t1 = -t16 * t26 + t25 * t64;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t41 - t38 * rSges(2,2) + r_base(1)) + g(2) * (t38 * rSges(2,1) + rSges(2,2) * t41 + r_base(2)) + g(3) * (rSges(2,3) + t57)) - m(3) * (g(1) * (rSges(3,1) * t17 - rSges(3,2) * t16 + t53) + g(2) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t56) + g(3) * (rSges(3,3) * t34 + t55) + (rSges(3,3) * t66 + g(3) * (rSges(3,1) * t37 + rSges(3,2) * t40) + (-rSges(3,3) - pkin(8)) * t65) * t32) - m(4) * (g(1) * (-rSges(4,2) * t17 + rSges(4,3) * t16 + t44) + g(2) * (-t15 * rSges(4,2) + t14 * rSges(4,3) + t45) + g(3) * (rSges(4,1) * t34 + t51) + (rSges(4,1) * t66 + g(3) * (-rSges(4,2) * t37 + (-rSges(4,3) - qJ(3)) * t40) + (-rSges(4,1) - pkin(8)) * t65) * t32) - m(5) * ((t47 * t66 + (-pkin(8) - t47) * t65) * t32 + (t51 + t47 * t34 + (t68 * t37 - t46 * t40) * t32) * g(3) + (t46 * t14 + t68 * t15 + t54) * g(2) + (t46 * t16 + t68 * t17 + t50) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t17 + t42) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t15 * rSges(6,3) + t43) + g(3) * (rSges(6,1) * t7 - rSges(6,2) * t6 + t48) + (g(3) * (t49 + (rSges(6,3) - t35) * t37) + t52) * t32) - m(7) * (g(1) * (t2 * pkin(5) + (t17 * t36 + t2 * t39) * rSges(7,1) + (t17 * t39 - t2 * t36) * rSges(7,2) + t69 * t1 + t42) + g(2) * (t4 * pkin(5) + (t15 * t36 + t39 * t4) * rSges(7,1) + (t15 * t39 - t36 * t4) * rSges(7,2) + t43 - t69 * t3) + t52 * t32 + (t69 * t6 + t48 + (t39 * rSges(7,1) - t36 * rSges(7,2) + pkin(5)) * t7 + (t49 + (t36 * rSges(7,1) + t39 * rSges(7,2) - t35) * t37) * t32) * g(3));
U  = t5;
