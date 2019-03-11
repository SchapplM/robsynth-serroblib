% Calculate Gravitation load on the joints for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
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
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:03
% EndTime: 2019-03-08 18:42:06
% DurationCPUTime: 0.81s
% Computational Cost: add. (771->127), mult. (2103->220), div. (0->0), fcn. (2705->16), ass. (0->74)
t43 = sin(pkin(12));
t44 = sin(pkin(11));
t48 = cos(pkin(12));
t49 = cos(pkin(11));
t51 = cos(pkin(6));
t81 = t49 * t51;
t32 = -t44 * t43 + t48 * t81;
t46 = sin(pkin(6));
t45 = sin(pkin(7));
t57 = cos(qJ(3));
t84 = t45 * t57;
t74 = t46 * t84;
t50 = cos(pkin(7));
t80 = t50 * t57;
t75 = pkin(3) * t80;
t33 = t43 * t81 + t44 * t48;
t54 = sin(qJ(3));
t90 = t33 * t54;
t58 = t32 * t75 + (-t49 * t74 - t90) * pkin(3);
t42 = sin(pkin(13));
t47 = cos(pkin(13));
t36 = t54 * t42 - t57 * t47;
t27 = t36 * t45;
t29 = t36 * t50;
t66 = t57 * t42 + t54 * t47;
t83 = t46 * t49;
t8 = t27 * t83 - t32 * t29 - t33 * t66;
t97 = t8 * pkin(4) + t58;
t96 = -t43 * t54 + t48 * t80;
t91 = rSges(7,3) + pkin(10);
t95 = g(1) * t91;
t28 = t66 * t45;
t30 = t66 * t50;
t19 = t51 * t28 + (t30 * t48 - t36 * t43) * t46;
t56 = cos(qJ(5));
t93 = t56 * pkin(5);
t92 = -rSges(6,3) - pkin(9);
t86 = t44 * t51;
t35 = -t43 * t86 + t49 * t48;
t89 = t35 * t54;
t87 = t44 * t46;
t85 = t45 * t46;
t82 = t46 * t50;
t52 = sin(qJ(6));
t78 = t52 * t56;
t55 = cos(qJ(6));
t77 = t55 * t56;
t76 = -m(5) - m(6) - m(7);
t72 = g(2) * t91;
t71 = g(3) * t91;
t69 = -m(3) - m(4) + t76;
t53 = sin(qJ(5));
t68 = rSges(6,1) * t56 - rSges(6,2) * t53;
t34 = -t49 * t43 - t48 * t86;
t65 = t34 * t75 + (t44 * t74 - t89) * pkin(3);
t63 = -t32 * t50 + t45 * t83;
t62 = t34 * t50 + t44 * t85;
t11 = -t27 * t87 - t34 * t29 - t35 * t66;
t61 = t11 * pkin(4) + t65;
t60 = (t96 * t46 + t51 * t84) * pkin(3);
t18 = -t51 * t27 + (-t29 * t48 - t43 * t66) * t46;
t59 = t18 * pkin(4) + t60;
t7 = t28 * t83 - t32 * t30 + t33 * t36;
t12 = t28 * t87 + t34 * t30 - t35 * t36;
t31 = -t48 * t85 + t51 * t50;
t21 = -t34 * t45 + t44 * t82;
t20 = -t32 * t45 - t49 * t82;
t14 = t19 * t56 + t31 * t53;
t13 = -t19 * t53 + t31 * t56;
t4 = t12 * t56 + t21 * t53;
t3 = -t12 * t53 + t21 * t56;
t2 = t20 * t53 - t56 * t7;
t1 = t20 * t56 + t53 * t7;
t5 = [(-m(2) + t69) * g(3), t69 * (g(3) * t51 + (g(1) * t44 - g(2) * t49) * t46) -m(4) * (g(1) * ((t62 * t57 - t89) * rSges(4,1) + (-t35 * t57 - t62 * t54) * rSges(4,2)) + g(2) * ((-t63 * t57 - t90) * rSges(4,1) + (-t33 * t57 + t63 * t54) * rSges(4,2)) + g(3) * ((rSges(4,1) * t57 - rSges(4,2) * t54) * t51 * t45 + (t96 * rSges(4,1) + (-t48 * t50 * t54 - t43 * t57) * rSges(4,2)) * t46)) - m(5) * (g(1) * (t11 * rSges(5,1) - rSges(5,2) * t12 + t65) + g(2) * (t8 * rSges(5,1) + t7 * rSges(5,2) + t58) + g(3) * (t18 * rSges(5,1) - rSges(5,2) * t19 + t60)) - m(6) * (g(1) * (t68 * t11 - t12 * t92 + t61) + g(2) * (t68 * t8 + t92 * t7 + t97) + g(3) * (t68 * t18 - t19 * t92 + t59)) - m(7) * (g(1) * (t11 * t93 + t12 * pkin(9) + (t11 * t77 + t12 * t52) * rSges(7,1) + (-t11 * t78 + t12 * t55) * rSges(7,2) + t61) + g(2) * (t8 * t93 - t7 * pkin(9) + (-t7 * t52 + t8 * t77) * rSges(7,1) + (-t7 * t55 - t8 * t78) * rSges(7,2) + t97) + g(3) * (t18 * t93 + t19 * pkin(9) + (t18 * t77 + t19 * t52) * rSges(7,1) + (-t18 * t78 + t19 * t55) * rSges(7,2) + t59) + (t11 * t95 + t18 * t71 + t8 * t72) * t53) t76 * (g(1) * t21 + g(2) * t20 + g(3) * t31) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (t13 * rSges(6,1) - t14 * rSges(6,2))) - m(7) * (t4 * t95 + t14 * t71 + t2 * t72 + (g(1) * t3 + g(2) * t1 + g(3) * t13) * (rSges(7,1) * t55 - rSges(7,2) * t52 + pkin(5))) -m(7) * (g(1) * ((-t11 * t55 - t4 * t52) * rSges(7,1) + (t11 * t52 - t4 * t55) * rSges(7,2)) + g(2) * ((-t2 * t52 - t8 * t55) * rSges(7,1) + (-t2 * t55 + t8 * t52) * rSges(7,2)) + g(3) * ((-t14 * t52 - t18 * t55) * rSges(7,1) + (-t14 * t55 + t18 * t52) * rSges(7,2)))];
taug  = t5(:);
