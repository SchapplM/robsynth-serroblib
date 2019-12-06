% Calculate potential energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:09
% EndTime: 2019-12-05 17:23:10
% DurationCPUTime: 0.45s
% Computational Cost: add. (374->116), mult. (889->159), div. (0->0), fcn. (1108->14), ass. (0->56)
t41 = sin(pkin(6));
t44 = cos(pkin(6));
t45 = cos(pkin(5));
t42 = sin(pkin(5));
t52 = cos(qJ(2));
t71 = t42 * t52;
t25 = -t41 * t71 + t45 * t44;
t40 = sin(pkin(11));
t43 = cos(pkin(11));
t49 = sin(qJ(2));
t68 = t45 * t52;
t28 = -t40 * t68 - t43 * t49;
t73 = t42 * t44;
t19 = -t28 * t41 + t40 * t73;
t79 = rSges(5,3) + pkin(9);
t78 = pkin(10) + rSges(6,3);
t77 = cos(qJ(3));
t75 = t40 * t42;
t74 = t42 * t43;
t72 = t42 * t49;
t69 = t45 * t49;
t67 = t40 * pkin(1) + r_base(2);
t64 = qJ(1) + r_base(3);
t63 = t41 * t77;
t62 = t44 * t77;
t61 = t43 * pkin(1) + pkin(7) * t75 + r_base(1);
t60 = t45 * pkin(7) + t64;
t59 = t42 * t63;
t26 = -t40 * t49 + t43 * t68;
t18 = -t26 * t41 - t43 * t73;
t29 = -t40 * t69 + t43 * t52;
t58 = t29 * pkin(2) + t19 * pkin(8) + t61;
t48 = sin(qJ(3));
t10 = t29 * t77 + (t28 * t44 + t41 * t75) * t48;
t57 = t10 * pkin(3) + t58;
t56 = pkin(2) * t72 + t25 * pkin(8) + t60;
t17 = t45 * t41 * t48 + (t44 * t48 * t52 + t49 * t77) * t42;
t55 = t17 * pkin(3) + t56;
t27 = t40 * t52 + t43 * t69;
t54 = t27 * pkin(2) - pkin(7) * t74 + pkin(8) * t18 + t67;
t8 = t27 * t77 + (t26 * t44 - t41 * t74) * t48;
t53 = t8 * pkin(3) + t54;
t51 = cos(qJ(4));
t50 = cos(qJ(5));
t47 = sin(qJ(4));
t46 = sin(qJ(5));
t16 = -t45 * t63 + t48 * t72 - t62 * t71;
t12 = t17 * t51 + t25 * t47;
t11 = t17 * t47 - t25 * t51;
t9 = -t28 * t62 + t29 * t48 - t40 * t59;
t7 = -t26 * t62 + t27 * t48 + t43 * t59;
t4 = t10 * t51 + t19 * t47;
t3 = t10 * t47 - t19 * t51;
t2 = t18 * t47 + t8 * t51;
t1 = -t18 * t51 + t8 * t47;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t43 * rSges(2,1) - t40 * rSges(2,2) + r_base(1)) + g(2) * (t40 * rSges(2,1) + t43 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t64)) - m(3) * (g(1) * (t29 * rSges(3,1) + t28 * rSges(3,2) + t61) + g(2) * (t27 * rSges(3,1) + t26 * rSges(3,2) + t67) + g(3) * (t45 * rSges(3,3) + t60) + (g(1) * rSges(3,3) * t40 + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t52) + g(2) * (-rSges(3,3) - pkin(7)) * t43) * t42) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t19 * rSges(4,3) + t58) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t18 * rSges(4,3) + t54) + g(3) * (t17 * rSges(4,1) - t16 * rSges(4,2) + t25 * rSges(4,3) + t56)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t79 * t9 + t57) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t7 * t79 + t53) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t79 * t16 + t55)) - m(6) * (g(1) * (t4 * pkin(4) + t9 * pkin(9) + (t4 * t50 + t9 * t46) * rSges(6,1) + (-t4 * t46 + t9 * t50) * rSges(6,2) + t78 * t3 + t57) + g(2) * (t2 * pkin(4) + t7 * pkin(9) + (t2 * t50 + t7 * t46) * rSges(6,1) + (-t2 * t46 + t7 * t50) * rSges(6,2) + t78 * t1 + t53) + g(3) * (t12 * pkin(4) + t16 * pkin(9) + (t12 * t50 + t16 * t46) * rSges(6,1) + (-t12 * t46 + t16 * t50) * rSges(6,2) + t78 * t11 + t55));
U = t5;
