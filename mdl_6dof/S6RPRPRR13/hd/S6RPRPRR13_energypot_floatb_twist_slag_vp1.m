% Calculate potential energy for
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:35
% EndTime: 2019-03-09 04:21:36
% DurationCPUTime: 0.49s
% Computational Cost: add. (498->135), mult. (1189->175), div. (0->0), fcn. (1483->14), ass. (0->64)
t49 = cos(pkin(6));
t44 = sin(pkin(12));
t56 = cos(qJ(1));
t78 = t56 * t44;
t47 = cos(pkin(12));
t53 = sin(qJ(1));
t79 = t53 * t47;
t31 = -t49 * t79 - t78;
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t46 = sin(pkin(6));
t84 = t46 * t53;
t22 = -t31 * t45 + t48 * t84;
t85 = t46 * t47;
t28 = -t45 * t85 + t49 * t48;
t91 = rSges(6,3) + pkin(10);
t90 = pkin(11) + rSges(7,3);
t89 = cos(qJ(3));
t87 = t44 * t46;
t52 = sin(qJ(3));
t86 = t45 * t52;
t83 = t46 * t56;
t82 = t48 * t52;
t80 = t53 * t44;
t77 = t56 * t47;
t76 = qJ(2) * t46;
t75 = rSges(5,3) + qJ(4);
t74 = pkin(8) + r_base(3);
t73 = t53 * pkin(1) + r_base(2);
t70 = t45 * t89;
t69 = t48 * t89;
t68 = t49 * qJ(2) + t74;
t67 = t56 * pkin(1) + t53 * t76 + r_base(1);
t66 = t46 * t70;
t29 = t49 * t77 - t80;
t21 = -t29 * t45 - t48 * t83;
t32 = -t49 * t80 + t77;
t65 = t32 * pkin(2) + t22 * pkin(9) + t67;
t64 = pkin(2) * t87 + t28 * pkin(9) + t68;
t14 = t32 * t89 + (t31 * t48 + t45 * t84) * t52;
t63 = t14 * pkin(3) + t65;
t18 = t49 * t86 + (t89 * t44 + t47 * t82) * t46;
t62 = t18 * pkin(3) + t64;
t13 = -t31 * t69 + t32 * t52 - t53 * t66;
t61 = t22 * pkin(4) + t13 * qJ(4) + t63;
t17 = -t49 * t70 + t52 * t87 - t69 * t85;
t60 = t28 * pkin(4) + t17 * qJ(4) + t62;
t30 = t49 * t78 + t79;
t59 = t30 * pkin(2) + t21 * pkin(9) - t56 * t76 + t73;
t12 = t29 * t82 + t30 * t89 - t83 * t86;
t58 = t12 * pkin(3) + t59;
t11 = -t29 * t69 + t30 * t52 + t56 * t66;
t57 = t21 * pkin(4) + t11 * qJ(4) + t58;
t55 = cos(qJ(5));
t54 = cos(qJ(6));
t51 = sin(qJ(5));
t50 = sin(qJ(6));
t10 = t17 * t51 + t28 * t55;
t9 = -t17 * t55 + t28 * t51;
t4 = t13 * t51 + t22 * t55;
t3 = -t13 * t55 + t22 * t51;
t2 = t11 * t51 + t21 * t55;
t1 = -t11 * t55 + t21 * t51;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t56 * rSges(2,1) - t53 * rSges(2,2) + r_base(1)) + g(2) * (t53 * rSges(2,1) + t56 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t74)) - m(3) * (g(1) * (t32 * rSges(3,1) + t31 * rSges(3,2) + t67) + g(2) * (t30 * rSges(3,1) + t29 * rSges(3,2) + t73) + g(3) * (t49 * rSges(3,3) + t68) + (g(1) * rSges(3,3) * t53 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t47) + g(2) * (-rSges(3,3) - qJ(2)) * t56) * t46) - m(4) * (g(1) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t22 * rSges(4,3) + t65) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t21 * rSges(4,3) + t59) + g(3) * (t18 * rSges(4,1) - t17 * rSges(4,2) + t28 * rSges(4,3) + t64)) - m(5) * (g(1) * (t22 * rSges(5,1) - t14 * rSges(5,2) + t75 * t13 + t63) + g(2) * (t21 * rSges(5,1) - t12 * rSges(5,2) + t75 * t11 + t58) + g(3) * (t28 * rSges(5,1) - t18 * rSges(5,2) + t75 * t17 + t62)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t91 * t14 + t61) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t91 * t12 + t57) + g(3) * (t10 * rSges(6,1) - t9 * rSges(6,2) + t91 * t18 + t60)) - m(7) * (g(1) * (t4 * pkin(5) + t14 * pkin(10) + (t14 * t50 + t4 * t54) * rSges(7,1) + (t14 * t54 - t4 * t50) * rSges(7,2) + t90 * t3 + t61) + g(2) * (t2 * pkin(5) + t12 * pkin(10) + (t12 * t50 + t2 * t54) * rSges(7,1) + (t12 * t54 - t2 * t50) * rSges(7,2) + t90 * t1 + t57) + g(3) * (t10 * pkin(5) + t18 * pkin(10) + (t10 * t54 + t18 * t50) * rSges(7,1) + (-t10 * t50 + t18 * t54) * rSges(7,2) + t90 * t9 + t60));
U  = t5;
