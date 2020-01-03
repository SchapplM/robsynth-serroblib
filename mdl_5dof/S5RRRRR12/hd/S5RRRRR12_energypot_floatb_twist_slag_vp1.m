% Calculate potential energy for
% S5RRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:29
% EndTime: 2019-12-31 22:46:29
% DurationCPUTime: 0.43s
% Computational Cost: add. (374->116), mult. (889->156), div. (0->0), fcn. (1108->14), ass. (0->57)
t43 = cos(pkin(5));
t47 = sin(qJ(2));
t52 = cos(qJ(1));
t69 = t52 * t47;
t48 = sin(qJ(1));
t51 = cos(qJ(2));
t70 = t48 * t51;
t28 = -t43 * t70 - t69;
t40 = sin(pkin(6));
t42 = cos(pkin(6));
t41 = sin(pkin(5));
t75 = t41 * t48;
t19 = -t28 * t40 + t42 * t75;
t74 = t41 * t51;
t25 = -t40 * t74 + t43 * t42;
t80 = rSges(5,3) + pkin(10);
t79 = pkin(11) + rSges(6,3);
t78 = cos(qJ(3));
t76 = t41 * t47;
t73 = t41 * t52;
t71 = t48 * t47;
t68 = t52 * t51;
t67 = pkin(7) + r_base(3);
t66 = t48 * pkin(1) + r_base(2);
t63 = t40 * t78;
t62 = t42 * t78;
t61 = t43 * pkin(8) + t67;
t60 = t52 * pkin(1) + pkin(8) * t75 + r_base(1);
t59 = t41 * t63;
t26 = t43 * t68 - t71;
t18 = -t26 * t40 - t42 * t73;
t29 = -t43 * t71 + t68;
t58 = t29 * pkin(2) + t19 * pkin(9) + t60;
t57 = pkin(2) * t76 + t25 * pkin(9) + t61;
t46 = sin(qJ(3));
t12 = t29 * t78 + (t28 * t42 + t40 * t75) * t46;
t56 = t12 * pkin(3) + t58;
t17 = t43 * t40 * t46 + (t42 * t46 * t51 + t78 * t47) * t41;
t55 = t17 * pkin(3) + t57;
t27 = t43 * t69 + t70;
t54 = t27 * pkin(2) - pkin(8) * t73 + t18 * pkin(9) + t66;
t10 = t27 * t78 + (t26 * t42 - t40 * t73) * t46;
t53 = t10 * pkin(3) + t54;
t50 = cos(qJ(4));
t49 = cos(qJ(5));
t45 = sin(qJ(4));
t44 = sin(qJ(5));
t16 = -t43 * t63 + t46 * t76 - t62 * t74;
t11 = -t28 * t62 + t29 * t46 - t48 * t59;
t9 = -t26 * t62 + t27 * t46 + t52 * t59;
t8 = t17 * t50 + t25 * t45;
t7 = t17 * t45 - t25 * t50;
t4 = t12 * t50 + t19 * t45;
t3 = t12 * t45 - t19 * t50;
t2 = t10 * t50 + t18 * t45;
t1 = t10 * t45 - t18 * t50;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t52 * rSges(2,1) - t48 * rSges(2,2) + r_base(1)) + g(2) * (t48 * rSges(2,1) + t52 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t67)) - m(3) * (g(1) * (t29 * rSges(3,1) + t28 * rSges(3,2) + t60) + g(2) * (t27 * rSges(3,1) + t26 * rSges(3,2) + t66) + g(3) * (t43 * rSges(3,3) + t61) + (g(1) * rSges(3,3) * t48 + g(3) * (rSges(3,1) * t47 + rSges(3,2) * t51) + g(2) * (-rSges(3,3) - pkin(8)) * t52) * t41) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t19 * rSges(4,3) + t58) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t18 * rSges(4,3) + t54) + g(3) * (t17 * rSges(4,1) - t16 * rSges(4,2) + t25 * rSges(4,3) + t57)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t80 * t11 + t56) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t80 * t9 + t53) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t80 * t16 + t55)) - m(6) * (g(1) * (t4 * pkin(4) + t11 * pkin(10) + (t11 * t44 + t4 * t49) * rSges(6,1) + (t11 * t49 - t4 * t44) * rSges(6,2) + t79 * t3 + t56) + g(2) * (t2 * pkin(4) + t9 * pkin(10) + (t2 * t49 + t9 * t44) * rSges(6,1) + (-t2 * t44 + t9 * t49) * rSges(6,2) + t79 * t1 + t53) + g(3) * (t8 * pkin(4) + t16 * pkin(10) + (t16 * t44 + t8 * t49) * rSges(6,1) + (t16 * t49 - t8 * t44) * rSges(6,2) + t79 * t7 + t55));
U = t5;
