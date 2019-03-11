% Calculate potential energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:07
% EndTime: 2019-03-09 04:01:08
% DurationCPUTime: 0.82s
% Computational Cost: add. (619->155), mult. (1464->210), div. (0->0), fcn. (1852->16), ass. (0->72)
t49 = sin(pkin(7));
t53 = cos(pkin(7));
t57 = sin(qJ(3));
t61 = cos(qJ(3));
t65 = t49 * (rSges(4,1) * t57 + rSges(4,2) * t61) + t53 * pkin(9);
t58 = sin(qJ(1));
t93 = g(1) * t58;
t62 = cos(qJ(1));
t92 = g(2) * t62;
t90 = -rSges(6,3) - pkin(10);
t89 = pkin(11) + rSges(7,3);
t52 = cos(pkin(12));
t48 = sin(pkin(12));
t80 = t58 * t48;
t54 = cos(pkin(6));
t81 = t54 * t62;
t35 = t52 * t81 - t80;
t88 = t35 * t49;
t79 = t58 * t52;
t37 = -t48 * t62 - t54 * t79;
t87 = t37 * t49;
t50 = sin(pkin(6));
t86 = t50 * t58;
t85 = t50 * t62;
t84 = t52 * t49;
t83 = t53 * t57;
t82 = t53 * t61;
t78 = pkin(9) + qJ(4);
t77 = pkin(8) + r_base(3);
t76 = t58 * pkin(1) + r_base(2);
t75 = t54 * qJ(2) + t77;
t74 = t62 * pkin(1) + qJ(2) * t86 + r_base(1);
t47 = sin(pkin(13));
t51 = cos(pkin(13));
t72 = t47 * t61 + t51 * t57;
t40 = -t47 * t57 + t51 * t61;
t32 = pkin(3) * t49 * t57 + t53 * t78;
t33 = pkin(3) * t83 - t49 * t78;
t43 = pkin(3) * t61 + pkin(2);
t71 = t54 * t32 + t75 + (t33 * t52 + t43 * t48) * t50;
t38 = t52 * t62 - t54 * t80;
t70 = t32 * t86 + t37 * t33 + t38 * t43 + t74;
t69 = t40 * t49;
t29 = t72 * t49;
t31 = t72 * t53;
t15 = t29 * t54 + (t31 * t52 + t40 * t48) * t50;
t68 = t15 * pkin(4) + t71;
t12 = t29 * t86 + t31 * t37 + t38 * t40;
t67 = t12 * pkin(4) + t70;
t66 = t50 * t69;
t36 = t48 * t81 + t79;
t64 = t35 * t33 + t36 * t43 + (-qJ(2) - t32) * t85 + t76;
t10 = -t29 * t85 + t35 * t31 + t36 * t40;
t63 = t10 * pkin(4) + t64;
t60 = cos(qJ(5));
t59 = cos(qJ(6));
t56 = sin(qJ(5));
t55 = sin(qJ(6));
t34 = -t50 * t84 + t53 * t54;
t30 = t40 * t53;
t21 = t53 * t86 - t87;
t20 = -t53 * t85 - t88;
t14 = (t30 * t52 - t48 * t72) * t50 + t54 * t69;
t11 = t37 * t30 - t38 * t72 + t58 * t66;
t9 = t30 * t35 - t36 * t72 - t62 * t66;
t8 = t15 * t60 + t34 * t56;
t7 = t15 * t56 - t34 * t60;
t4 = t12 * t60 + t21 * t56;
t3 = t12 * t56 - t21 * t60;
t2 = t10 * t60 + t20 * t56;
t1 = t10 * t56 - t20 * t60;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t62 - t58 * rSges(2,2) + r_base(1)) + g(2) * (t58 * rSges(2,1) + rSges(2,2) * t62 + r_base(2)) + g(3) * (rSges(2,3) + t77)) - m(3) * (g(1) * (rSges(3,1) * t38 + rSges(3,2) * t37 + t74) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + t76) + g(3) * (rSges(3,3) * t54 + t75) + (rSges(3,3) * t93 + g(3) * (rSges(3,1) * t48 + rSges(3,2) * t52) + (-rSges(3,3) - qJ(2)) * t92) * t50) - m(4) * (g(1) * (t38 * pkin(2) - pkin(9) * t87 + (t37 * t83 + t38 * t61) * rSges(4,1) + (t37 * t82 - t38 * t57) * rSges(4,2) + t21 * rSges(4,3) + t74) + g(2) * (t36 * pkin(2) - pkin(9) * t88 + (t35 * t83 + t36 * t61) * rSges(4,1) + (t35 * t82 - t36 * t57) * rSges(4,2) + t20 * rSges(4,3) + t76) + (t65 * t93 + (-qJ(2) - t65) * t92) * t50 + (t34 * rSges(4,3) + t75 + t65 * t54 + (t48 * pkin(2) - pkin(9) * t84 + (t48 * t61 + t52 * t83) * rSges(4,1) + (-t48 * t57 + t52 * t82) * rSges(4,2)) * t50) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t12 + rSges(5,2) * t11 + rSges(5,3) * t21 + t70) + g(2) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t20 * rSges(5,3) + t64) + g(3) * (rSges(5,1) * t15 + rSges(5,2) * t14 + rSges(5,3) * t34 + t71)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t90 * t11 + t67) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t90 * t9 + t63) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t90 * t14 + t68)) - m(7) * (g(1) * (t4 * pkin(5) - t11 * pkin(10) + (-t11 * t55 + t4 * t59) * rSges(7,1) + (-t11 * t59 - t4 * t55) * rSges(7,2) + t89 * t3 + t67) + g(2) * (t2 * pkin(5) - t9 * pkin(10) + (t2 * t59 - t55 * t9) * rSges(7,1) + (-t2 * t55 - t59 * t9) * rSges(7,2) + t89 * t1 + t63) + g(3) * (t8 * pkin(5) - t14 * pkin(10) + (-t14 * t55 + t59 * t8) * rSges(7,1) + (-t14 * t59 - t55 * t8) * rSges(7,2) + t89 * t7 + t68));
U  = t5;
