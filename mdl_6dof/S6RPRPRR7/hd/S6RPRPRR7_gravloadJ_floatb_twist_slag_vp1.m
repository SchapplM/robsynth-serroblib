% Calculate Gravitation load on the joints for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:58
% EndTime: 2019-03-09 03:54:59
% DurationCPUTime: 0.54s
% Computational Cost: add. (322->109), mult. (323->141), div. (0->0), fcn. (277->10), ass. (0->57)
t32 = qJ(3) + pkin(10);
t27 = qJ(5) + t32;
t24 = cos(t27);
t77 = rSges(7,3) + pkin(9);
t78 = t77 * t24;
t39 = cos(qJ(1));
t70 = g(2) * t39;
t36 = sin(qJ(1));
t76 = g(1) * t36 - t70;
t75 = g(1) * t39 + g(2) * t36;
t34 = sin(qJ(6));
t64 = rSges(7,2) * t34;
t51 = t24 * t64;
t74 = t51 * t70;
t35 = sin(qJ(3));
t69 = t35 * pkin(3);
t38 = cos(qJ(3));
t68 = t38 * pkin(3);
t67 = rSges(4,3) + pkin(7);
t65 = rSges(6,1) * t24;
t23 = sin(t27);
t63 = t23 * t36;
t62 = t23 * t39;
t61 = t24 * t36;
t60 = t36 * t34;
t37 = cos(qJ(6));
t59 = t36 * t37;
t58 = t39 * t34;
t57 = t39 * t37;
t33 = -qJ(4) - pkin(7);
t56 = rSges(5,3) - t33;
t55 = pkin(1) * t39 + qJ(2) * t36;
t54 = -m(5) - m(6) - m(7);
t25 = sin(t32);
t10 = pkin(4) * t25 + t69;
t29 = t39 * qJ(2);
t31 = -pkin(8) + t33;
t53 = t10 * t39 + t31 * t36 + t29;
t52 = t10 * t36 + t55;
t50 = -rSges(7,1) * t37 - pkin(5);
t49 = rSges(6,1) * t61 - rSges(6,2) * t63;
t47 = rSges(4,1) * t35 + rSges(4,2) * t38;
t46 = rSges(6,1) * t23 + rSges(6,2) * t24;
t45 = g(1) * t29 + g(2) * t55;
t26 = cos(t32);
t44 = rSges(5,1) * t25 + rSges(5,2) * t26 + t69;
t42 = t24 * rSges(7,1) * t59 + pkin(5) * t61 - t36 * t51 + t63 * t77;
t41 = t78 + (t50 + t64) * t23;
t40 = -t23 * t77 + t24 * t50;
t15 = rSges(6,2) * t62;
t11 = pkin(4) * t26 + t68;
t6 = t36 * t11;
t4 = t23 * t57 - t60;
t3 = t23 * t58 + t59;
t2 = t23 * t59 + t58;
t1 = -t23 * t60 + t57;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t36 - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 - rSges(2,2) * t36)) - m(3) * (g(1) * (t39 * rSges(3,3) + t29 + (rSges(3,2) - pkin(1)) * t36) + g(2) * (-rSges(3,2) * t39 + rSges(3,3) * t36 + t55)) - m(4) * ((g(1) * t47 + g(2) * t67) * t39 + (g(1) * (-pkin(1) - t67) + g(2) * t47) * t36 + t45) - m(5) * ((g(1) * t44 + g(2) * t56) * t39 + (g(1) * (-pkin(1) - t56) + g(2) * t44) * t36 + t45) - m(6) * (g(1) * t53 + g(2) * t52 + (g(1) * t46 + g(2) * (rSges(6,3) - t31)) * t39 + (g(1) * (-rSges(6,3) - pkin(1)) + g(2) * t46) * t36) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) - t36 * pkin(1) + pkin(5) * t62 + t53) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t63 - t39 * t31 + t52) - t75 * t78) (-m(3) - m(4) + t54) * t76, -m(4) * (-g(3) * t47 + t76 * (rSges(4,1) * t38 - rSges(4,2) * t35)) - m(5) * (-g(3) * t44 + t76 * (rSges(5,1) * t26 - rSges(5,2) * t25 + t68)) - m(6) * (g(1) * (t49 + t6) + g(2) * (t15 + (-t11 - t65) * t39) + g(3) * (-t10 - t46)) - m(7) * (g(1) * (t42 + t6) + t74 + g(3) * (-t10 + t41) + (-t11 + t40) * t70) t54 * t75, -m(6) * (g(1) * t49 + g(2) * (-t39 * t65 + t15) - g(3) * t46) - m(7) * (g(1) * t42 + g(3) * t41 + t40 * t70 + t74) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t4) + g(3) * (-rSges(7,1) * t34 - rSges(7,2) * t37) * t24)];
taug  = t5(:);
