% Calculate potential energy for
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:23
% EndTime: 2019-03-09 15:14:23
% DurationCPUTime: 0.34s
% Computational Cost: add. (313->120), mult. (668->148), div. (0->0), fcn. (786->10), ass. (0->54)
t76 = rSges(7,1) + pkin(5);
t38 = sin(pkin(6));
t44 = cos(qJ(2));
t75 = t38 * t44;
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t74 = t40 * t41;
t42 = sin(qJ(1));
t73 = t41 * t42;
t45 = cos(qJ(1));
t72 = t41 * t45;
t71 = t42 * t44;
t70 = t45 * t40;
t43 = cos(qJ(3));
t69 = t45 * t43;
t68 = qJ(4) * t38;
t39 = cos(pkin(6));
t67 = qJ(4) * t39;
t66 = rSges(7,2) + qJ(5);
t65 = rSges(6,3) + qJ(5);
t64 = rSges(7,3) + qJ(6);
t63 = cos(pkin(10));
t62 = pkin(7) + r_base(3);
t61 = t42 * pkin(1) + r_base(2);
t60 = t38 * t74;
t59 = t41 * t67;
t58 = t41 * pkin(2) + t62;
t57 = t39 * t63;
t56 = t41 * t63;
t55 = t45 * pkin(1) + t42 * pkin(8) + r_base(1);
t54 = t38 * t56;
t53 = t45 * t44 * pkin(2) + pkin(9) * t72 + t55;
t52 = pkin(2) * t71 - t45 * pkin(8) + pkin(9) * t73 + t61;
t51 = qJ(4) * t60 + t41 * t43 * pkin(3) + (-pkin(9) - t67) * t44 + t58;
t19 = t42 * t43 - t44 * t70;
t20 = t42 * t40 + t44 * t69;
t50 = t20 * pkin(3) - t19 * t68 + t45 * t59 + t53;
t37 = sin(pkin(10));
t9 = t43 * t56 + (-t39 * t74 - t75) * t37;
t49 = t9 * pkin(4) + t51;
t6 = t20 * t63 + (t19 * t39 + t38 * t72) * t37;
t48 = t6 * pkin(4) + t50;
t17 = -t40 * t71 - t69;
t18 = t43 * t71 - t70;
t47 = t18 * pkin(3) - t17 * t68 + t42 * t59 + t52;
t4 = t18 * t63 + (t17 * t39 + t38 * t73) * t37;
t46 = t4 * pkin(4) + t47;
t16 = -t44 * t39 + t60;
t11 = -t19 * t38 + t39 * t72;
t10 = -t17 * t38 + t39 * t73;
t8 = t63 * t75 + (t37 * t43 + t40 * t57) * t41;
t5 = -t19 * t57 + t20 * t37 - t45 * t54;
t3 = -t17 * t57 + t18 * t37 - t42 * t54;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t45 * rSges(2,1) - t42 * rSges(2,2) + r_base(1)) + g(2) * (t42 * rSges(2,1) + t45 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t62)) - m(3) * (g(1) * (t42 * rSges(3,3) + t55) + g(2) * (rSges(3,1) * t71 - rSges(3,2) * t73 + t61) + g(3) * (t41 * rSges(3,1) + t44 * rSges(3,2) + t62) + (g(1) * (rSges(3,1) * t44 - rSges(3,2) * t41) + g(2) * (-rSges(3,3) - pkin(8))) * t45) - m(4) * (g(1) * (t20 * rSges(4,1) + t19 * rSges(4,2) + rSges(4,3) * t72 + t53) + g(2) * (t18 * rSges(4,1) + t17 * rSges(4,2) + rSges(4,3) * t73 + t52) + g(3) * ((-rSges(4,3) - pkin(9)) * t44 + (rSges(4,1) * t43 - rSges(4,2) * t40) * t41 + t58)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t11 * rSges(5,3) + t50) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t10 * rSges(5,3) + t47) + g(3) * (t9 * rSges(5,1) - t8 * rSges(5,2) + t16 * rSges(5,3) + t51)) - m(6) * (g(1) * (t11 * rSges(6,1) - t6 * rSges(6,2) + t65 * t5 + t48) + g(2) * (t10 * rSges(6,1) - t4 * rSges(6,2) + t65 * t3 + t46) + g(3) * (t16 * rSges(6,1) - t9 * rSges(6,2) + t65 * t8 + t49)) - m(7) * (g(1) * (t76 * t11 + t66 * t5 + t64 * t6 + t48) + g(2) * (t76 * t10 + t66 * t3 + t64 * t4 + t46) + g(3) * (t76 * t16 + t64 * t9 + t66 * t8 + t49));
U  = t1;
