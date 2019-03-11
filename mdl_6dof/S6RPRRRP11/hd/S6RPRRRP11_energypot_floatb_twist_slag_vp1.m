% Calculate potential energy for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:21
% EndTime: 2019-03-09 06:34:21
% DurationCPUTime: 0.51s
% Computational Cost: add. (572->137), mult. (1388->178), div. (0->0), fcn. (1753->14), ass. (0->67)
t49 = cos(pkin(6));
t44 = sin(pkin(12));
t56 = cos(qJ(1));
t79 = t56 * t44;
t47 = cos(pkin(12));
t54 = sin(qJ(1));
t80 = t54 * t47;
t31 = -t49 * t80 - t79;
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t46 = sin(pkin(6));
t84 = t46 * t54;
t64 = -t31 * t45 + t48 * t84;
t85 = t46 * t47;
t63 = -t45 * t85 + t49 * t48;
t91 = rSges(5,3) + pkin(10);
t90 = rSges(6,3) + pkin(11);
t89 = cos(qJ(3));
t88 = cos(qJ(4));
t86 = t44 * t46;
t83 = t46 * t56;
t81 = t54 * t44;
t78 = t56 * t47;
t77 = rSges(7,3) + qJ(6) + pkin(11);
t76 = t46 * qJ(2);
t75 = pkin(8) + r_base(3);
t74 = t54 * pkin(1) + r_base(2);
t51 = sin(qJ(5));
t71 = pkin(5) * t51 + pkin(10);
t70 = t45 * t89;
t69 = t48 * t89;
t68 = t49 * qJ(2) + t75;
t67 = t56 * pkin(1) + t54 * t76 + r_base(1);
t66 = t46 * t70;
t29 = t49 * t78 - t81;
t65 = -t29 * t45 - t48 * t83;
t32 = -t49 * t81 + t78;
t62 = t32 * pkin(2) + t64 * pkin(9) + t67;
t61 = pkin(2) * t86 + t63 * pkin(9) + t68;
t53 = sin(qJ(3));
t18 = t32 * t89 + (t31 * t48 + t45 * t84) * t53;
t60 = t18 * pkin(3) + t62;
t23 = t49 * t45 * t53 + (t47 * t48 * t53 + t89 * t44) * t46;
t59 = t23 * pkin(3) + t61;
t30 = t49 * t79 + t80;
t58 = t30 * pkin(2) + t65 * pkin(9) - t56 * t76 + t74;
t16 = t30 * t89 + (t29 * t48 - t45 * t83) * t53;
t57 = t16 * pkin(3) + t58;
t55 = cos(qJ(5));
t52 = sin(qJ(4));
t40 = t55 * pkin(5) + pkin(4);
t22 = -t49 * t70 + t53 * t86 - t69 * t85;
t17 = -t31 * t69 + t32 * t53 - t54 * t66;
t15 = -t29 * t69 + t30 * t53 + t56 * t66;
t14 = t23 * t88 + t63 * t52;
t13 = t23 * t52 - t63 * t88;
t10 = t18 * t88 + t64 * t52;
t9 = t18 * t52 - t64 * t88;
t8 = t16 * t88 + t65 * t52;
t7 = t16 * t52 - t65 * t88;
t6 = t14 * t55 + t22 * t51;
t5 = -t14 * t51 + t22 * t55;
t4 = t10 * t55 + t17 * t51;
t3 = -t10 * t51 + t17 * t55;
t2 = t15 * t51 + t8 * t55;
t1 = t15 * t55 - t8 * t51;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t56 * rSges(2,1) - t54 * rSges(2,2) + r_base(1)) + g(2) * (t54 * rSges(2,1) + t56 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t75)) - m(3) * (g(1) * (t32 * rSges(3,1) + t31 * rSges(3,2) + t67) + g(2) * (t30 * rSges(3,1) + t29 * rSges(3,2) + t74) + g(3) * (t49 * rSges(3,3) + t68) + (g(1) * rSges(3,3) * t54 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t47) + g(2) * (-rSges(3,3) - qJ(2)) * t56) * t46) - m(4) * (g(1) * (t18 * rSges(4,1) - t17 * rSges(4,2) + t64 * rSges(4,3) + t62) + g(2) * (t16 * rSges(4,1) - t15 * rSges(4,2) + t65 * rSges(4,3) + t58) + g(3) * (t23 * rSges(4,1) - t22 * rSges(4,2) + t63 * rSges(4,3) + t61)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t91 * t17 + t60) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t91 * t15 + t57) + g(3) * (t14 * rSges(5,1) - t13 * rSges(5,2) + t91 * t22 + t59)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * pkin(4) + t17 * pkin(10) + t90 * t9 + t60) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) + t15 * pkin(10) + t90 * t7 + t57) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t14 * pkin(4) + t22 * pkin(10) + t90 * t13 + t59)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t10 * t40 + t71 * t17 + t77 * t9 + t60) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t71 * t15 + t8 * t40 + t77 * t7 + t57) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t77 * t13 + t14 * t40 + t71 * t22 + t59));
U  = t11;
