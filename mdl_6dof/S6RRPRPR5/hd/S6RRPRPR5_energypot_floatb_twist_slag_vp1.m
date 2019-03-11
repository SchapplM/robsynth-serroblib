% Calculate potential energy for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:51
% EndTime: 2019-03-09 10:28:52
% DurationCPUTime: 0.55s
% Computational Cost: add. (421->129), mult. (853->159), div. (0->0), fcn. (1059->14), ass. (0->54)
t36 = sin(pkin(11));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t59 = cos(pkin(11));
t21 = -t42 * t36 + t45 * t59;
t71 = pkin(2) * t42;
t69 = rSges(5,3) + pkin(9);
t37 = sin(pkin(6));
t43 = sin(qJ(1));
t67 = t43 * t37;
t66 = t43 * t42;
t65 = t43 * t45;
t46 = cos(qJ(1));
t64 = t46 * t37;
t63 = t46 * t42;
t62 = t46 * t45;
t61 = pkin(10) + qJ(5) + rSges(7,3);
t60 = qJ(5) + rSges(6,3);
t58 = pkin(7) + r_base(3);
t29 = t45 * pkin(2) + pkin(1);
t57 = t46 * t29 + r_base(1);
t39 = cos(pkin(6));
t56 = t39 * pkin(8) + t58;
t19 = t39 * t71 + (-pkin(8) - qJ(3)) * t37;
t54 = t46 * t19 + t43 * t29 + r_base(2);
t20 = -t45 * t36 - t42 * t59;
t18 = t20 * t39;
t8 = -t46 * t18 + t43 * t21;
t53 = t8 * pkin(3) + t54;
t52 = t39 * qJ(3) + t37 * t71 + t56;
t34 = pkin(12) + qJ(6);
t30 = sin(t34);
t31 = cos(t34);
t38 = cos(pkin(12));
t51 = t31 * rSges(7,1) - t30 * rSges(7,2) + t38 * pkin(5) + pkin(4);
t10 = t43 * t18 + t46 * t21;
t50 = t10 * pkin(3) - t43 * t19 + t57;
t17 = t20 * t37;
t49 = -t17 * pkin(3) + t52;
t35 = sin(pkin(12));
t48 = t30 * rSges(7,1) + t31 * rSges(7,2) + t35 * pkin(5) + pkin(9);
t47 = t21 * t39;
t44 = cos(qJ(4));
t41 = sin(qJ(4));
t16 = t21 * t37;
t12 = -t17 * t44 + t39 * t41;
t11 = -t17 * t41 - t39 * t44;
t9 = t46 * t20 - t43 * t47;
t7 = t43 * t20 + t46 * t47;
t4 = t10 * t44 + t41 * t67;
t3 = t10 * t41 - t44 * t67;
t2 = -t41 * t64 + t8 * t44;
t1 = t8 * t41 + t44 * t64;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t46 * rSges(2,1) - t43 * rSges(2,2) + r_base(1)) + g(2) * (t43 * rSges(2,1) + t46 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t58)) - m(3) * (g(1) * (t46 * pkin(1) + r_base(1) + (-t39 * t66 + t62) * rSges(3,1) + (-t39 * t65 - t63) * rSges(3,2)) + g(2) * (t43 * pkin(1) + r_base(2) + (t39 * t63 + t65) * rSges(3,1) + (t39 * t62 - t66) * rSges(3,2)) + g(3) * (t39 * rSges(3,3) + t56) + (g(3) * (rSges(3,1) * t42 + rSges(3,2) * t45) + (g(1) * t43 - g(2) * t46) * (rSges(3,3) + pkin(8))) * t37) - m(4) * (g(1) * (t10 * rSges(4,1) + t9 * rSges(4,2) + (rSges(4,3) * t37 - t19) * t43 + t57) + g(2) * (t8 * rSges(4,1) + t7 * rSges(4,2) - rSges(4,3) * t64 + t54) + g(3) * (-t17 * rSges(4,1) + t16 * rSges(4,2) + t39 * rSges(4,3) + t52)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t69 * t9 + t50) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) - t69 * t7 + t53) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) - t69 * t16 + t49)) - m(6) * (g(1) * (t4 * pkin(4) - t9 * pkin(9) + (-t9 * t35 + t4 * t38) * rSges(6,1) + (-t4 * t35 - t9 * t38) * rSges(6,2) + t60 * t3 + t50) + g(2) * (t2 * pkin(4) - t7 * pkin(9) + (t2 * t38 - t7 * t35) * rSges(6,1) + (-t2 * t35 - t7 * t38) * rSges(6,2) + t60 * t1 + t53) + g(3) * (t12 * pkin(4) - t16 * pkin(9) + (t12 * t38 - t16 * t35) * rSges(6,1) + (-t12 * t35 - t16 * t38) * rSges(6,2) + t60 * t11 + t49)) - m(7) * (g(1) * (t61 * t3 + t51 * t4 - t48 * t9 + t50) + g(2) * (t61 * t1 + t51 * t2 - t48 * t7 + t53) + g(3) * (t61 * t11 + t51 * t12 - t48 * t16 + t49));
U  = t5;
