% Calculate Gravitation load on the joints for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:35
% EndTime: 2019-03-09 08:33:37
% DurationCPUTime: 0.69s
% Computational Cost: add. (232->135), mult. (453->180), div. (0->0), fcn. (415->6), ass. (0->56)
t69 = rSges(7,1) + pkin(5);
t61 = rSges(6,3) + pkin(8);
t50 = rSges(7,3) + qJ(6) + pkin(8);
t27 = cos(qJ(1));
t24 = sin(qJ(1));
t64 = g(2) * t24;
t34 = g(1) * t27 + t64;
t68 = -pkin(2) - pkin(3);
t23 = sin(qJ(2));
t67 = pkin(4) * t23;
t22 = sin(qJ(5));
t66 = pkin(5) * t22;
t26 = cos(qJ(2));
t63 = g(3) * t26;
t18 = t26 * pkin(2);
t60 = rSges(4,1) * t26;
t59 = rSges(5,1) * t23;
t58 = rSges(5,2) * t26;
t25 = cos(qJ(5));
t14 = pkin(5) * t25 + pkin(4);
t57 = t14 * t23;
t56 = t23 * t27;
t55 = t24 * t22;
t54 = t24 * t25;
t53 = t24 * t26;
t52 = t26 * t27;
t51 = t27 * t25;
t15 = t23 * qJ(3);
t49 = t15 + t18;
t48 = t27 * pkin(1) + t24 * pkin(7);
t47 = qJ(3) * t26;
t46 = -rSges(5,3) - qJ(4);
t45 = -m(5) - m(6) - m(7);
t44 = rSges(5,2) + t68;
t43 = t26 * pkin(3) + t49;
t42 = -t61 + t68;
t41 = -t50 + t68;
t40 = -pkin(1) - t15;
t39 = -qJ(4) - t66;
t38 = pkin(2) * t52 + t27 * t15 + t48;
t37 = g(1) * t44;
t36 = pkin(3) * t52 + t38;
t35 = g(1) * t42;
t33 = g(1) * t41;
t32 = rSges(3,1) * t26 - rSges(3,2) * t23;
t30 = rSges(6,1) * t25 - rSges(6,2) * t22 + pkin(4);
t4 = -t22 * t56 - t54;
t2 = t23 * t55 - t51;
t29 = rSges(7,1) * t25 - rSges(7,2) * t22 + t14;
t7 = t24 * t47;
t9 = t27 * t47;
t28 = g(1) * t9 + g(2) * t7 + g(3) * t43;
t19 = t27 * pkin(7);
t5 = t23 * t51 - t55;
t3 = -t22 * t27 - t23 * t54;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t24 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t27 + t19) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t56 + t48) + (g(1) * (-pkin(1) - t32) + g(2) * rSges(3,3)) * t24) - m(4) * (g(1) * (rSges(4,2) * t27 + t19) + g(2) * (rSges(4,1) * t52 + rSges(4,3) * t56 + t38) + (g(1) * (-rSges(4,3) * t23 - t18 + t40 - t60) + g(2) * rSges(4,2)) * t24) - m(5) * (g(1) * t19 + g(2) * t36 + (g(1) * t46 + g(2) * (-t58 + t59)) * t27 + (g(1) * (t40 - t59) + g(2) * t46 + t26 * t37) * t24) - m(6) * (g(1) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t19) + g(2) * (t5 * rSges(6,1) + t4 * rSges(6,2) + t36) + (-g(1) * qJ(4) + g(2) * (t61 * t26 + t67)) * t27 + (g(1) * (t40 - t67) - g(2) * qJ(4) + t26 * t35) * t24) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t19) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t36) + (g(1) * t39 + g(2) * (t50 * t26 + t57)) * t27 + (g(1) * (t40 - t57) + g(2) * t39 + t26 * t33) * t24) -m(3) * (g(3) * t32 + t34 * (-rSges(3,1) * t23 - rSges(3,2) * t26)) - m(4) * (g(1) * (rSges(4,3) * t52 + t9) + g(2) * (rSges(4,3) * t53 + t7) + g(3) * (t49 + t60) + (g(3) * rSges(4,3) + t34 * (-rSges(4,1) - pkin(2))) * t23) - m(5) * (g(1) * (rSges(5,1) * t52 + t9) + g(2) * (rSges(5,1) * t53 + t7) + g(3) * (t43 - t58) + (g(3) * rSges(5,1) + t27 * t37 + t44 * t64) * t23) - m(6) * ((g(3) * t61 + t34 * t30) * t26 + (g(3) * t30 + t27 * t35 + t42 * t64) * t23 + t28) - m(7) * ((g(3) * t50 + t34 * t29) * t26 + (g(3) * t29 + t27 * t33 + t41 * t64) * t23 + t28) (-m(4) + t45) * (t23 * t34 - t63) t45 * (-g(1) * t24 + g(2) * t27) -m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 + rSges(6,2) * t3)) - m(7) * (g(1) * (-t5 * rSges(7,2) + t69 * t4) + g(2) * (t3 * rSges(7,2) - t69 * t2)) + (-m(6) * (rSges(6,1) * t22 + rSges(6,2) * t25) - m(7) * (rSges(7,1) * t22 + rSges(7,2) * t25 + t66)) * t63, -m(7) * (g(3) * t23 + t26 * t34)];
taug  = t1(:);
