% Calculate potential energy for
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:40
% EndTime: 2019-03-09 11:57:41
% DurationCPUTime: 0.55s
% Computational Cost: add. (409->130), mult. (853->161), div. (0->0), fcn. (1059->12), ass. (0->57)
t38 = sin(pkin(11));
t44 = sin(qJ(2));
t48 = cos(qJ(2));
t61 = cos(pkin(11));
t27 = -t44 * t38 + t48 * t61;
t73 = pkin(2) * t44;
t71 = rSges(5,3) + pkin(9);
t70 = rSges(6,3) + pkin(10);
t39 = sin(pkin(6));
t45 = sin(qJ(1));
t68 = t45 * t39;
t67 = t45 * t44;
t66 = t45 * t48;
t49 = cos(qJ(1));
t65 = t49 * t39;
t64 = t49 * t44;
t63 = t49 * t48;
t62 = rSges(7,3) + qJ(6) + pkin(10);
t60 = pkin(7) + r_base(3);
t35 = t48 * pkin(2) + pkin(1);
t59 = t49 * t35 + r_base(1);
t42 = sin(qJ(5));
t58 = pkin(5) * t42 + pkin(9);
t40 = cos(pkin(6));
t57 = t40 * pkin(8) + t60;
t25 = t40 * t73 + (-pkin(8) - qJ(3)) * t39;
t55 = t49 * t25 + t45 * t35 + r_base(2);
t26 = -t48 * t38 - t44 * t61;
t24 = t26 * t40;
t14 = -t49 * t24 + t45 * t27;
t54 = t14 * pkin(3) + t55;
t53 = t40 * qJ(3) + t39 * t73 + t57;
t16 = t45 * t24 + t49 * t27;
t52 = t16 * pkin(3) - t45 * t25 + t59;
t23 = t26 * t39;
t51 = -t23 * pkin(3) + t53;
t50 = t27 * t40;
t47 = cos(qJ(4));
t46 = cos(qJ(5));
t43 = sin(qJ(4));
t34 = t46 * pkin(5) + pkin(4);
t22 = t27 * t39;
t18 = -t23 * t47 + t40 * t43;
t17 = -t23 * t43 - t40 * t47;
t15 = t49 * t26 - t45 * t50;
t13 = t45 * t26 + t49 * t50;
t10 = t16 * t47 + t43 * t68;
t9 = t16 * t43 - t47 * t68;
t8 = t14 * t47 - t43 * t65;
t7 = t14 * t43 + t47 * t65;
t6 = t18 * t46 - t22 * t42;
t5 = -t18 * t42 - t22 * t46;
t4 = t10 * t46 - t15 * t42;
t3 = -t10 * t42 - t15 * t46;
t2 = -t13 * t42 + t8 * t46;
t1 = -t13 * t46 - t8 * t42;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t49 * rSges(2,1) - t45 * rSges(2,2) + r_base(1)) + g(2) * (t45 * rSges(2,1) + t49 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t60)) - m(3) * (g(1) * (t49 * pkin(1) + r_base(1) + (-t40 * t67 + t63) * rSges(3,1) + (-t40 * t66 - t64) * rSges(3,2)) + g(2) * (t45 * pkin(1) + r_base(2) + (t40 * t64 + t66) * rSges(3,1) + (t40 * t63 - t67) * rSges(3,2)) + g(3) * (t40 * rSges(3,3) + t57) + (g(3) * (rSges(3,1) * t44 + rSges(3,2) * t48) + (g(1) * t45 - g(2) * t49) * (rSges(3,3) + pkin(8))) * t39) - m(4) * (g(1) * (t16 * rSges(4,1) + t15 * rSges(4,2) + (rSges(4,3) * t39 - t25) * t45 + t59) + g(2) * (t14 * rSges(4,1) + t13 * rSges(4,2) - rSges(4,3) * t65 + t55) + g(3) * (-t23 * rSges(4,1) + t22 * rSges(4,2) + t40 * rSges(4,3) + t53)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) - t71 * t15 + t52) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) - t71 * t13 + t54) + g(3) * (t18 * rSges(5,1) - t17 * rSges(5,2) - t71 * t22 + t51)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * pkin(4) - t15 * pkin(9) + t70 * t9 + t52) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) - t13 * pkin(9) + t70 * t7 + t54) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t18 * pkin(4) - t22 * pkin(9) + t70 * t17 + t51)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t10 * t34 - t58 * t15 + t62 * t9 + t52) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t58 * t13 + t8 * t34 + t62 * t7 + t54) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t17 * t62 + t18 * t34 - t58 * t22 + t51));
U  = t11;
