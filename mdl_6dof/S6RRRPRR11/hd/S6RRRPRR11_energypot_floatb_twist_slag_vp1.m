% Calculate potential energy for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:43
% EndTime: 2019-03-09 19:24:43
% DurationCPUTime: 0.44s
% Computational Cost: add. (355->124), mult. (737->148), div. (0->0), fcn. (900->12), ass. (0->57)
t76 = rSges(5,2) + pkin(9);
t75 = rSges(4,3) + pkin(9);
t74 = rSges(6,3) - pkin(9);
t73 = pkin(11) + rSges(7,3);
t72 = cos(qJ(3));
t37 = sin(pkin(6));
t41 = sin(qJ(2));
t71 = t37 * t41;
t42 = sin(qJ(1));
t70 = t37 * t42;
t45 = cos(qJ(2));
t69 = t37 * t45;
t46 = cos(qJ(1));
t68 = t37 * t46;
t67 = rSges(5,3) + qJ(4);
t66 = cos(pkin(6));
t65 = pkin(7) + r_base(3);
t64 = -pkin(10) - t74;
t63 = t42 * pkin(1) + r_base(2);
t62 = t37 * t72;
t61 = t66 * pkin(8) + t65;
t60 = t42 * t66;
t59 = t46 * t66;
t58 = t46 * pkin(1) + pkin(8) * t70 + r_base(1);
t57 = pkin(2) * t71 + t61;
t27 = -t41 * t60 + t46 * t45;
t56 = t27 * pkin(2) + t58;
t40 = sin(qJ(3));
t23 = t66 * t40 + t41 * t62;
t55 = t23 * pkin(3) + t57;
t16 = t27 * t72 + t40 * t70;
t54 = t16 * pkin(3) + t56;
t38 = sin(qJ(6));
t43 = cos(qJ(6));
t53 = t43 * rSges(7,1) - t38 * rSges(7,2) + pkin(5);
t52 = -t38 * rSges(7,1) - t43 * rSges(7,2) + pkin(9) - pkin(10);
t25 = t41 * t59 + t42 * t45;
t51 = t25 * pkin(2) - pkin(8) * t68 + t63;
t14 = t25 * t72 - t40 * t68;
t50 = t14 * pkin(3) + t51;
t15 = t27 * t40 - t42 * t62;
t49 = t16 * pkin(4) + t15 * qJ(4) + t54;
t22 = t40 * t71 - t66 * t72;
t48 = t23 * pkin(4) + pkin(10) * t69 + t22 * qJ(4) + t55;
t13 = t25 * t40 + t46 * t62;
t47 = t14 * pkin(4) + t13 * qJ(4) + t50;
t44 = cos(qJ(5));
t39 = sin(qJ(5));
t26 = t46 * t41 + t45 * t60;
t24 = t42 * t41 - t45 * t59;
t6 = t22 * t39 + t23 * t44;
t5 = -t22 * t44 + t23 * t39;
t4 = t15 * t39 + t16 * t44;
t3 = -t15 * t44 + t16 * t39;
t2 = t13 * t39 + t14 * t44;
t1 = -t13 * t44 + t14 * t39;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t46 * rSges(2,1) - t42 * rSges(2,2) + r_base(1)) + g(2) * (t42 * rSges(2,1) + t46 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t65)) - m(3) * (g(1) * (t27 * rSges(3,1) - t26 * rSges(3,2) + t58) + g(2) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t63) + g(3) * (t66 * rSges(3,3) + t61) + (g(1) * rSges(3,3) * t42 + g(3) * (rSges(3,1) * t41 + rSges(3,2) * t45) + g(2) * (-rSges(3,3) - pkin(8)) * t46) * t37) - m(4) * (g(1) * (t16 * rSges(4,1) - t15 * rSges(4,2) + t75 * t26 + t56) + g(2) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t75 * t24 + t51) + g(3) * (t23 * rSges(4,1) - t22 * rSges(4,2) - t75 * t69 + t57)) - m(5) * (g(1) * (t16 * rSges(5,1) + t67 * t15 + t76 * t26 + t54) + g(2) * (t14 * rSges(5,1) + t67 * t13 + t76 * t24 + t50) + g(3) * (t23 * rSges(5,1) + t67 * t22 - t76 * t69 + t55)) - m(6) * (g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t74 * t69 + t48) + (t2 * rSges(6,1) - t1 * rSges(6,2) + t64 * t24 + t47) * g(2) + (t4 * rSges(6,1) - t3 * rSges(6,2) + t64 * t26 + t49) * g(1)) - m(7) * (g(1) * (t52 * t26 + t73 * t3 + t53 * t4 + t49) + g(2) * (t73 * t1 + t53 * t2 + t52 * t24 + t47) + g(3) * (t6 * pkin(5) - pkin(9) * t69 + (t38 * t69 + t6 * t43) * rSges(7,1) + (-t6 * t38 + t43 * t69) * rSges(7,2) + t73 * t5 + t48));
U  = t7;
