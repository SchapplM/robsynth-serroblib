% Calculate potential energy for
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR15_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR15_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:25
% EndTime: 2019-03-09 20:24:26
% DurationCPUTime: 0.51s
% Computational Cost: add. (498->135), mult. (1189->174), div. (0->0), fcn. (1483->14), ass. (0->63)
t47 = cos(pkin(6));
t51 = sin(qJ(2));
t56 = cos(qJ(1));
t77 = t56 * t51;
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t78 = t52 * t55;
t31 = -t47 * t78 - t77;
t44 = sin(pkin(7));
t46 = cos(pkin(7));
t45 = sin(pkin(6));
t84 = t45 * t52;
t22 = -t31 * t44 + t46 * t84;
t83 = t45 * t55;
t28 = -t44 * t83 + t47 * t46;
t90 = rSges(6,3) + pkin(11);
t89 = pkin(12) + rSges(7,3);
t88 = cos(qJ(3));
t50 = sin(qJ(3));
t86 = t44 * t50;
t85 = t45 * t51;
t82 = t45 * t56;
t81 = t46 * t50;
t79 = t52 * t51;
t76 = t56 * t55;
t75 = rSges(5,3) + qJ(4);
t74 = pkin(8) + r_base(3);
t73 = t52 * pkin(1) + r_base(2);
t70 = t44 * t88;
t69 = t46 * t88;
t68 = t47 * pkin(9) + t74;
t67 = t56 * pkin(1) + pkin(9) * t84 + r_base(1);
t66 = t45 * t70;
t29 = t47 * t76 - t79;
t21 = -t29 * t44 - t46 * t82;
t32 = -t47 * t79 + t76;
t65 = t32 * pkin(2) + t22 * pkin(10) + t67;
t64 = pkin(2) * t85 + t28 * pkin(10) + t68;
t14 = t32 * t88 + (t31 * t46 + t44 * t84) * t50;
t63 = t14 * pkin(3) + t65;
t18 = t47 * t86 + (t88 * t51 + t55 * t81) * t45;
t62 = t18 * pkin(3) + t64;
t13 = -t31 * t69 + t32 * t50 - t52 * t66;
t61 = t22 * pkin(4) + t13 * qJ(4) + t63;
t17 = -t47 * t70 + t50 * t85 - t69 * t83;
t60 = t28 * pkin(4) + t17 * qJ(4) + t62;
t30 = t47 * t77 + t78;
t59 = t30 * pkin(2) - pkin(9) * t82 + t21 * pkin(10) + t73;
t12 = t29 * t81 + t30 * t88 - t82 * t86;
t58 = t12 * pkin(3) + t59;
t11 = -t29 * t69 + t30 * t50 + t56 * t66;
t57 = t21 * pkin(4) + t11 * qJ(4) + t58;
t54 = cos(qJ(5));
t53 = cos(qJ(6));
t49 = sin(qJ(5));
t48 = sin(qJ(6));
t10 = t17 * t49 + t28 * t54;
t9 = -t17 * t54 + t28 * t49;
t4 = t13 * t49 + t22 * t54;
t3 = -t13 * t54 + t22 * t49;
t2 = t11 * t49 + t21 * t54;
t1 = -t11 * t54 + t21 * t49;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t56 * rSges(2,1) - t52 * rSges(2,2) + r_base(1)) + g(2) * (t52 * rSges(2,1) + t56 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t74)) - m(3) * (g(1) * (t32 * rSges(3,1) + t31 * rSges(3,2) + t67) + g(2) * (t30 * rSges(3,1) + t29 * rSges(3,2) + t73) + g(3) * (t47 * rSges(3,3) + t68) + (g(1) * rSges(3,3) * t52 + g(3) * (rSges(3,1) * t51 + rSges(3,2) * t55) + g(2) * (-rSges(3,3) - pkin(9)) * t56) * t45) - m(4) * (g(1) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t22 * rSges(4,3) + t65) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t21 * rSges(4,3) + t59) + g(3) * (t18 * rSges(4,1) - t17 * rSges(4,2) + t28 * rSges(4,3) + t64)) - m(5) * (g(1) * (t22 * rSges(5,1) - t14 * rSges(5,2) + t75 * t13 + t63) + g(2) * (t21 * rSges(5,1) - t12 * rSges(5,2) + t75 * t11 + t58) + g(3) * (t28 * rSges(5,1) - t18 * rSges(5,2) + t75 * t17 + t62)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t90 * t14 + t61) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t90 * t12 + t57) + g(3) * (t10 * rSges(6,1) - t9 * rSges(6,2) + t90 * t18 + t60)) - m(7) * (g(1) * (t4 * pkin(5) + t14 * pkin(11) + (t14 * t48 + t4 * t53) * rSges(7,1) + (t14 * t53 - t4 * t48) * rSges(7,2) + t89 * t3 + t61) + g(2) * (t2 * pkin(5) + t12 * pkin(11) + (t12 * t48 + t2 * t53) * rSges(7,1) + (t12 * t53 - t2 * t48) * rSges(7,2) + t89 * t1 + t57) + g(3) * (t10 * pkin(5) + t18 * pkin(11) + (t10 * t53 + t18 * t48) * rSges(7,1) + (-t10 * t48 + t18 * t53) * rSges(7,2) + t89 * t9 + t60));
U  = t5;
