% Calculate potential energy for
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:41
% EndTime: 2019-03-09 21:41:41
% DurationCPUTime: 0.44s
% Computational Cost: add. (356->116), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->53)
t72 = rSges(7,2) + qJ(5);
t71 = rSges(7,3) + qJ(6);
t70 = rSges(6,1) + pkin(10);
t69 = rSges(4,3) + pkin(9);
t68 = rSges(5,3) + pkin(10);
t67 = cos(qJ(3));
t36 = sin(pkin(6));
t39 = sin(qJ(2));
t66 = t36 * t39;
t40 = sin(qJ(1));
t65 = t36 * t40;
t42 = cos(qJ(2));
t64 = t36 * t42;
t43 = cos(qJ(1));
t63 = t36 * t43;
t62 = rSges(6,3) + qJ(5);
t61 = cos(pkin(6));
t60 = pkin(7) + r_base(3);
t58 = t40 * pkin(1) + r_base(2);
t57 = t36 * t67;
t56 = t61 * pkin(8) + t60;
t55 = t40 * t61;
t54 = t43 * t61;
t53 = t43 * pkin(1) + pkin(8) * t65 + r_base(1);
t52 = pkin(2) * t66 + t56;
t26 = -t39 * t55 + t43 * t42;
t51 = t26 * pkin(2) + t53;
t24 = t39 * t54 + t40 * t42;
t50 = t24 * pkin(2) - pkin(8) * t63 + t58;
t38 = sin(qJ(3));
t15 = t26 * t67 + t38 * t65;
t25 = t43 * t39 + t42 * t55;
t49 = t15 * pkin(3) + t25 * pkin(9) + t51;
t37 = sin(qJ(4));
t41 = cos(qJ(4));
t6 = t15 * t41 + t25 * t37;
t48 = t6 * pkin(4) + t49;
t22 = t38 * t61 + t39 * t57;
t47 = t22 * pkin(3) - pkin(9) * t64 + t52;
t11 = t22 * t41 - t37 * t64;
t46 = t11 * pkin(4) + t47;
t13 = t24 * t67 - t38 * t63;
t23 = t40 * t39 - t42 * t54;
t45 = t13 * pkin(3) + t23 * pkin(9) + t50;
t4 = t13 * t41 + t23 * t37;
t44 = t4 * pkin(4) + t45;
t21 = t38 * t66 - t61 * t67;
t14 = t26 * t38 - t40 * t57;
t12 = t24 * t38 + t43 * t57;
t10 = t22 * t37 + t41 * t64;
t5 = t15 * t37 - t25 * t41;
t3 = t13 * t37 - t23 * t41;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t43 * rSges(2,1) - t40 * rSges(2,2) + r_base(1)) + g(2) * (t40 * rSges(2,1) + t43 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t60)) - m(3) * (g(1) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t53) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t58) + g(3) * (t61 * rSges(3,3) + t56) + (g(1) * rSges(3,3) * t40 + g(3) * (rSges(3,1) * t39 + rSges(3,2) * t42) + g(2) * (-rSges(3,3) - pkin(8)) * t43) * t36) - m(4) * (g(1) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t25 * t69 + t51) + g(2) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t23 * t69 + t50) + g(3) * (t22 * rSges(4,1) - t21 * rSges(4,2) - t64 * t69 + t52)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t14 * t68 + t49) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t12 * t68 + t45) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t21 * t68 + t47)) - m(6) * (g(1) * (-t6 * rSges(6,2) + t14 * t70 + t62 * t5 + t48) + g(2) * (-t4 * rSges(6,2) + t12 * t70 + t62 * t3 + t44) + g(3) * (-t11 * rSges(6,2) + t62 * t10 + t21 * t70 + t46)) - m(7) * (g(1) * (t72 * t5 + t71 * t6 + t48) + g(2) * (t72 * t3 + t71 * t4 + t44) + g(3) * (t72 * t10 + t71 * t11 + t46) + (g(1) * t14 + g(2) * t12 + g(3) * t21) * (rSges(7,1) + pkin(5) + pkin(10)));
U  = t1;
