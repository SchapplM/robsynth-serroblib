% Calculate Gravitation load on the joints for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:52
% EndTime: 2019-03-10 03:27:53
% DurationCPUTime: 0.71s
% Computational Cost: add. (711->130), mult. (498->158), div. (0->0), fcn. (428->12), ass. (0->77)
t107 = rSges(7,3) + pkin(11);
t42 = qJ(2) + qJ(3);
t38 = qJ(4) + t42;
t35 = qJ(5) + t38;
t28 = sin(t35);
t29 = cos(t35);
t43 = sin(qJ(6));
t90 = rSges(7,2) * t43;
t106 = t29 * t107 + t28 * t90;
t104 = t29 * rSges(6,1) - rSges(6,2) * t28;
t32 = sin(t38);
t33 = cos(t38);
t75 = t33 * rSges(5,1) - rSges(5,2) * t32;
t36 = sin(t42);
t37 = cos(t42);
t76 = t37 * rSges(4,1) - rSges(4,2) * t36;
t68 = -rSges(5,1) * t32 - rSges(5,2) * t33;
t99 = pkin(3) * t36;
t103 = t68 - t99;
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t102 = g(1) * t48 + g(2) * t45;
t58 = t29 * pkin(5) + t107 * t28;
t46 = cos(qJ(6));
t94 = rSges(7,1) * t46;
t101 = (-pkin(5) - t94) * t28;
t49 = -pkin(8) - pkin(7);
t44 = sin(qJ(2));
t100 = pkin(2) * t44;
t98 = pkin(4) * t32;
t95 = rSges(3,3) + pkin(7);
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t34 = t40 + pkin(1);
t87 = t43 * t45;
t86 = t43 * t48;
t85 = t45 * t46;
t84 = t46 * t48;
t83 = rSges(4,3) - t49;
t41 = -pkin(9) + t49;
t82 = rSges(5,3) - t41;
t39 = -pkin(10) + t41;
t81 = rSges(6,3) - t39;
t31 = pkin(3) * t37;
t15 = t31 + t34;
t79 = t106 * t45;
t78 = t106 * t48;
t73 = t31 + t75;
t27 = pkin(4) * t33;
t72 = t27 + t104;
t12 = -t98 - t99;
t71 = rSges(3,1) * t47 - rSges(3,2) * t44;
t69 = -rSges(4,1) * t36 - rSges(4,2) * t37;
t66 = -rSges(6,1) * t28 - rSges(6,2) * t29;
t65 = t31 + t72;
t64 = pkin(1) + t71;
t63 = t34 + t76;
t62 = t15 + t75;
t60 = t66 * t45;
t59 = t66 * t48;
t57 = t58 + (-t90 + t94) * t29;
t54 = t27 + t57;
t53 = t31 + t54;
t52 = g(1) * t78 + g(2) * t79;
t51 = t102 * t101;
t11 = t48 * t12;
t10 = t45 * t12;
t9 = t12 - t100;
t8 = t27 + t15;
t7 = t29 * t84 + t87;
t6 = -t29 * t86 + t85;
t5 = -t29 * t85 + t86;
t4 = t29 * t87 + t84;
t3 = t48 * t9;
t2 = t45 * t9;
t1 = t48 * t8;
t13 = [-m(2) * (g(1) * (-rSges(2,1) * t45 - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - rSges(2,2) * t45)) - m(3) * ((g(1) * t95 + g(2) * t64) * t48 + (-g(1) * t64 + g(2) * t95) * t45) - m(4) * ((g(1) * t83 + g(2) * t63) * t48 + (-g(1) * t63 + g(2) * t83) * t45) - m(5) * ((g(1) * t82 + g(2) * t62) * t48 + (-g(1) * t62 + g(2) * t82) * t45) - m(6) * (g(2) * t1 + (g(1) * t81 + g(2) * t104) * t48 + (g(1) * (-t104 - t8) + g(2) * t81) * t45) - m(7) * (g(1) * (rSges(7,1) * t5 + rSges(7,2) * t4) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 + t1) + (-g(1) * t39 + g(2) * t58) * t48 + (g(1) * (-t58 - t8) - g(2) * t39) * t45) -m(3) * (g(3) * t71 + t102 * (-rSges(3,1) * t44 - rSges(3,2) * t47)) - m(4) * (g(3) * (t40 + t76) + t102 * (t69 - t100)) - m(5) * (g(3) * (t40 + t73) + t102 * (-t100 + t103)) - m(6) * (g(1) * (t3 + t59) + g(2) * (t2 + t60) + g(3) * (t40 + t65)) - m(7) * (g(1) * (t3 + t78) + g(2) * (t2 + t79) + g(3) * (t40 + t53) + t51) -m(4) * (g(3) * t76 + t102 * t69) - m(5) * (g(3) * t73 + t102 * t103) - m(6) * (g(1) * (t11 + t59) + g(2) * (t10 + t60) + g(3) * t65) - m(7) * (g(1) * (t11 + t78) + g(2) * (t10 + t79) + g(3) * t53 + t51) -m(7) * t52 + (-m(5) * t75 - m(6) * t72 - m(7) * t54) * g(3) + t102 * (-m(5) * t68 - m(6) * (t66 - t98) - m(7) * (-t98 + t101)) -m(6) * (g(1) * t59 + g(2) * t60 + g(3) * t104) - m(7) * (g(3) * t57 + t51 + t52) -m(7) * (g(1) * (rSges(7,1) * t6 - rSges(7,2) * t7) + g(2) * (-rSges(7,1) * t4 + rSges(7,2) * t5) + g(3) * (-rSges(7,1) * t43 - rSges(7,2) * t46) * t28)];
taug  = t13(:);
