% Calculate potential energy for
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:41:54
% EndTime: 2019-03-08 18:41:55
% DurationCPUTime: 0.82s
% Computational Cost: add. (619->155), mult. (1464->212), div. (0->0), fcn. (1852->16), ass. (0->72)
t50 = sin(pkin(7));
t55 = cos(pkin(7));
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t65 = (rSges(4,1) * t59 + rSges(4,2) * t62) * t50 + t55 * pkin(8);
t49 = sin(pkin(11));
t93 = g(1) * t49;
t54 = cos(pkin(11));
t92 = g(2) * t54;
t90 = -rSges(6,3) - pkin(9);
t89 = pkin(10) + rSges(7,3);
t48 = sin(pkin(12));
t53 = cos(pkin(12));
t56 = cos(pkin(6));
t81 = t54 * t56;
t35 = -t48 * t49 + t53 * t81;
t88 = t35 * t50;
t85 = t49 * t56;
t37 = -t48 * t54 - t53 * t85;
t87 = t37 * t50;
t51 = sin(pkin(6));
t86 = t49 * t51;
t84 = t51 * t54;
t83 = t51 * t55;
t82 = t53 * t50;
t80 = t55 * t59;
t79 = t55 * t62;
t78 = pkin(8) + qJ(4);
t77 = t49 * pkin(1) + r_base(2);
t76 = qJ(1) + r_base(3);
t75 = t54 * pkin(1) + qJ(2) * t86 + r_base(1);
t74 = t56 * qJ(2) + t76;
t47 = sin(pkin(13));
t52 = cos(pkin(13));
t72 = t47 * t62 + t59 * t52;
t40 = -t59 * t47 + t52 * t62;
t32 = pkin(3) * t50 * t59 + t78 * t55;
t33 = pkin(3) * t80 - t78 * t50;
t38 = -t48 * t85 + t53 * t54;
t43 = pkin(3) * t62 + pkin(2);
t71 = t32 * t86 + t37 * t33 + t38 * t43 + t75;
t70 = t40 * t50;
t69 = t56 * t32 + t74 + (t33 * t53 + t43 * t48) * t51;
t29 = t72 * t50;
t31 = t72 * t55;
t10 = t29 * t86 + t31 * t37 + t38 * t40;
t68 = t10 * pkin(4) + t71;
t67 = t51 * t70;
t15 = t29 * t56 + (t31 * t53 + t40 * t48) * t51;
t66 = t15 * pkin(4) + t69;
t36 = t48 * t81 + t49 * t53;
t64 = t35 * t33 + t36 * t43 + (-qJ(2) - t32) * t84 + t77;
t8 = -t29 * t84 + t31 * t35 + t36 * t40;
t63 = t8 * pkin(4) + t64;
t61 = cos(qJ(5));
t60 = cos(qJ(6));
t58 = sin(qJ(5));
t57 = sin(qJ(6));
t34 = -t51 * t82 + t55 * t56;
t30 = t40 * t55;
t21 = t49 * t83 - t87;
t20 = -t54 * t83 - t88;
t14 = (t30 * t53 - t48 * t72) * t51 + t56 * t70;
t12 = t15 * t61 + t34 * t58;
t11 = t15 * t58 - t34 * t61;
t9 = t37 * t30 - t38 * t72 + t49 * t67;
t7 = t30 * t35 - t36 * t72 - t54 * t67;
t4 = t10 * t61 + t21 * t58;
t3 = t10 * t58 - t21 * t61;
t2 = t20 * t58 + t61 * t8;
t1 = -t20 * t61 + t58 * t8;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t54 - rSges(2,2) * t49 + r_base(1)) + g(2) * (rSges(2,1) * t49 + rSges(2,2) * t54 + r_base(2)) + g(3) * (rSges(2,3) + t76)) - m(3) * (g(1) * (rSges(3,1) * t38 + rSges(3,2) * t37 + t75) + g(2) * (rSges(3,1) * t36 + rSges(3,2) * t35 + t77) + g(3) * (rSges(3,3) * t56 + t74) + (rSges(3,3) * t93 + g(3) * (rSges(3,1) * t48 + rSges(3,2) * t53) + (-rSges(3,3) - qJ(2)) * t92) * t51) - m(4) * (g(1) * (t38 * pkin(2) - pkin(8) * t87 + (t37 * t80 + t38 * t62) * rSges(4,1) + (t37 * t79 - t38 * t59) * rSges(4,2) + t21 * rSges(4,3) + t75) + g(2) * (t36 * pkin(2) - pkin(8) * t88 + (t35 * t80 + t36 * t62) * rSges(4,1) + (t35 * t79 - t36 * t59) * rSges(4,2) + t20 * rSges(4,3) + t77) + (t65 * t93 + (-qJ(2) - t65) * t92) * t51 + (t34 * rSges(4,3) + t74 + t65 * t56 + (t48 * pkin(2) - pkin(8) * t82 + (t48 * t62 + t53 * t80) * rSges(4,1) + (-t48 * t59 + t53 * t79) * rSges(4,2)) * t51) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t10 + rSges(5,2) * t9 + rSges(5,3) * t21 + t71) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t7 + rSges(5,3) * t20 + t64) + g(3) * (rSges(5,1) * t15 + rSges(5,2) * t14 + rSges(5,3) * t34 + t69)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t90 * t9 + t68) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t90 * t7 + t63) + g(3) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t90 * t14 + t66)) - m(7) * (g(1) * (t4 * pkin(5) - t9 * pkin(9) + (t4 * t60 - t57 * t9) * rSges(7,1) + (-t4 * t57 - t60 * t9) * rSges(7,2) + t89 * t3 + t68) + g(2) * (t2 * pkin(5) - t7 * pkin(9) + (t2 * t60 - t57 * t7) * rSges(7,1) + (-t2 * t57 - t60 * t7) * rSges(7,2) + t89 * t1 + t63) + g(3) * (t12 * pkin(5) - t14 * pkin(9) + (t12 * t60 - t14 * t57) * rSges(7,1) + (-t12 * t57 - t14 * t60) * rSges(7,2) + t89 * t11 + t66));
U  = t5;
