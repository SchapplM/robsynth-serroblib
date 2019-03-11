% Calculate potential energy for
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:07
% EndTime: 2019-03-09 09:07:07
% DurationCPUTime: 0.54s
% Computational Cost: add. (253->124), mult. (469->142), div. (0->0), fcn. (530->10), ass. (0->54)
t29 = cos(pkin(6));
t33 = sin(qJ(1));
t36 = cos(qJ(2));
t58 = t33 * t36;
t32 = sin(qJ(2));
t37 = cos(qJ(1));
t59 = t32 * t37;
t14 = t29 * t59 + t58;
t28 = sin(pkin(6));
t55 = qJ(4) * t28;
t69 = t14 * pkin(3) + t37 * t55;
t68 = g(1) * t33;
t67 = g(2) * t37;
t66 = pkin(10) + rSges(7,3);
t65 = t28 * t32;
t64 = t28 * t33;
t35 = cos(qJ(5));
t63 = t28 * t35;
t62 = t28 * t36;
t61 = t28 * t37;
t60 = t32 * t33;
t57 = t36 * t37;
t56 = qJ(3) * t36;
t54 = rSges(6,3) - qJ(3);
t53 = pkin(7) + r_base(3);
t52 = t33 * pkin(1) + r_base(2);
t51 = -pkin(9) - t54;
t50 = t29 * pkin(8) + t53;
t49 = t14 * pkin(2) + t52;
t48 = t37 * pkin(1) + pkin(8) * t64 + r_base(1);
t47 = pkin(2) * t65 + t50;
t16 = -t29 * t60 + t57;
t46 = t16 * pkin(2) + t48;
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t45 = t34 * rSges(7,1) - t30 * rSges(7,2) + pkin(5);
t13 = -t29 * t57 + t60;
t44 = t13 * qJ(3) + t49;
t43 = -rSges(7,1) * t30 - rSges(7,2) * t34 - pkin(9) + qJ(3);
t15 = t29 * t58 + t59;
t42 = t15 * qJ(3) + t46;
t41 = pkin(3) * t65 - t29 * qJ(4) + t47;
t40 = t14 * pkin(4) - pkin(8) * t61 + t49 + t69;
t9 = t16 * pkin(3);
t39 = t16 * pkin(4) - t33 * t55 + t46 + t9;
t38 = pkin(4) * t65 + pkin(9) * t62 + t41;
t31 = sin(qJ(5));
t12 = -t29 * t31 + t32 * t63;
t11 = t29 * t35 + t31 * t65;
t4 = t16 * t35 - t31 * t64;
t3 = t16 * t31 + t33 * t63;
t2 = t14 * t35 + t31 * t61;
t1 = t14 * t31 - t35 * t61;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t37 - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + rSges(2,2) * t37 + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t48) + g(2) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t52) + g(3) * (rSges(3,3) * t29 + t50) + (rSges(3,3) * t68 + g(3) * (rSges(3,1) * t32 + rSges(3,2) * t36) + (-rSges(3,3) - pkin(8)) * t67) * t28) - m(4) * (g(1) * (rSges(4,1) * t16 + rSges(4,3) * t15 + t42) + g(2) * (t14 * rSges(4,1) + t13 * rSges(4,3) + t44) + g(3) * (rSges(4,2) * t29 + t47) + (rSges(4,2) * t68 + g(3) * (rSges(4,1) * t32 - rSges(4,3) * t36 - t56) + (-rSges(4,2) - pkin(8)) * t67) * t28) - m(5) * (g(1) * (rSges(5,1) * t16 + rSges(5,2) * t15 + t42 + t9) + g(2) * (t14 * rSges(5,1) + t13 * rSges(5,2) + t44 + t69) + g(3) * (-rSges(5,3) * t29 + t41) + (g(3) * (rSges(5,1) * t32 + (-rSges(5,2) - qJ(3)) * t36) + (rSges(5,3) - pkin(8)) * t67 + (-rSges(5,3) - qJ(4)) * t68) * t28) - m(6) * (g(3) * (t12 * rSges(6,1) - rSges(6,2) * t11 + t54 * t62 + t38) + (t2 * rSges(6,1) - t1 * rSges(6,2) + t51 * t13 + t40) * g(2) + (rSges(6,1) * t4 - rSges(6,2) * t3 + t51 * t15 + t39) * g(1)) - m(7) * (g(1) * (t43 * t15 + t66 * t3 + t45 * t4 + t39) + g(2) * (t66 * t1 + t43 * t13 + t45 * t2 + t40) + g(3) * (t12 * pkin(5) - t28 * t56 + (t12 * t34 + t30 * t62) * rSges(7,1) + (-t12 * t30 + t34 * t62) * rSges(7,2) + t66 * t11 + t38));
U  = t5;
