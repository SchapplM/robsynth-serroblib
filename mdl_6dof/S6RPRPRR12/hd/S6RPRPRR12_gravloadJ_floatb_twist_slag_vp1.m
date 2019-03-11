% Calculate Gravitation load on the joints for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:13
% EndTime: 2019-03-09 04:18:15
% DurationCPUTime: 0.62s
% Computational Cost: add. (237->124), mult. (401->166), div. (0->0), fcn. (367->8), ass. (0->55)
t30 = sin(qJ(3));
t67 = g(3) * t30;
t65 = rSges(6,3) + pkin(8);
t55 = rSges(7,3) + pkin(9) + pkin(8);
t33 = cos(qJ(3));
t24 = t33 * qJ(4);
t72 = -t33 * rSges(5,3) - t24;
t34 = cos(qJ(1));
t68 = g(2) * t34;
t31 = sin(qJ(1));
t69 = g(1) * t31;
t46 = -t68 + t69;
t71 = -pkin(1) - pkin(7);
t29 = sin(qJ(5));
t70 = pkin(5) * t29;
t66 = rSges(5,2) - pkin(3);
t64 = t30 * t31;
t63 = t30 * t34;
t62 = t31 * t33;
t61 = t33 * rSges(4,2);
t28 = qJ(5) + qJ(6);
t21 = sin(t28);
t59 = t34 * t21;
t22 = cos(t28);
t58 = t34 * t22;
t57 = t34 * t29;
t32 = cos(qJ(5));
t56 = t34 * t32;
t54 = pkin(3) * t62 + qJ(4) * t64;
t25 = t34 * qJ(2);
t53 = pkin(3) * t63 + t25;
t52 = t34 * pkin(1) + t31 * qJ(2);
t51 = -m(5) - m(6) - m(7);
t50 = -pkin(3) - t65;
t49 = -pkin(3) - t55;
t48 = t34 * pkin(7) + t52;
t47 = pkin(3) * t64 + t48;
t44 = t30 * rSges(4,1) + t61;
t43 = rSges(6,1) * t29 + rSges(6,2) * t32;
t42 = -t31 * t29 + t33 * t56;
t12 = -t32 * t62 - t57;
t41 = g(1) * t54 + g(3) * t24;
t40 = rSges(7,1) * t21 + rSges(7,2) * t22 + t70;
t39 = t65 * t30 - t24;
t38 = -t30 * rSges(5,2) + t72;
t37 = t55 * t30 - t33 * t70 - t24;
t5 = t31 * t21 - t33 * t58;
t6 = -t31 * t22 - t33 * t59;
t7 = -t22 * t62 - t59;
t8 = -t21 * t62 + t58;
t36 = g(1) * (t7 * rSges(7,1) - t8 * rSges(7,2)) + g(2) * (-t5 * rSges(7,1) + t6 * rSges(7,2)) + (rSges(7,1) * t22 - rSges(7,2) * t21) * t67;
t20 = t32 * pkin(5) + pkin(4);
t13 = -t29 * t62 + t56;
t11 = -t31 * t32 - t33 * t57;
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) - t31 * rSges(2,2))) - m(3) * (g(1) * (t34 * rSges(3,3) + t25 + (rSges(3,2) - pkin(1)) * t31) + g(2) * (-t34 * rSges(3,2) + t31 * rSges(3,3) + t52)) - m(4) * (g(1) * (rSges(4,1) * t63 + t34 * t61 + t25) + g(2) * (t34 * rSges(4,3) + t48) + (g(1) * (-rSges(4,3) + t71) + g(2) * t44) * t31) - m(5) * (g(1) * t53 + g(2) * t47 + (g(2) * rSges(5,1) + g(1) * t38) * t34 + (g(1) * (-rSges(5,1) + t71) + g(2) * t38) * t31) - m(6) * (g(1) * (t11 * rSges(6,1) - rSges(6,2) * t42 + t53) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t47) + (g(2) * pkin(4) + g(1) * t39) * t34 + (g(1) * (-pkin(4) + t71) + g(2) * t39) * t31) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t53) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t47) + (g(1) * t37 + g(2) * t20) * t34 + (g(1) * (-t20 + t71) + g(2) * t37) * t31) (-m(3) - m(4) + t51) * t46, -m(4) * (-g(3) * t44 + t46 * (rSges(4,1) * t33 - rSges(4,2) * t30)) - m(5) * (g(1) * ((-rSges(5,2) * t33 + rSges(5,3) * t30) * t31 + t54) + g(3) * (t66 * t30 - t72) + (t66 * t33 + (-rSges(5,3) - qJ(4)) * t30) * t68) - m(6) * ((g(3) * t43 + t65 * t69) * t33 + (g(3) * t50 + t43 * t69) * t30 + (t50 * t33 + (-qJ(4) - t43) * t30) * t68 + t41) - m(7) * ((g(3) * t40 + t49 * t68 + t55 * t69) * t33 + (g(3) * t49 + t40 * t69 + (-qJ(4) - t40) * t68) * t30 + t41) t51 * (-t46 * t33 + t67) -m(6) * (g(1) * (t12 * rSges(6,1) - t13 * rSges(6,2)) + g(2) * (rSges(6,1) * t42 + t11 * rSges(6,2)) + (rSges(6,1) * t32 - rSges(6,2) * t29) * t67) - m(7) * ((g(1) * t12 + g(2) * t42 + t32 * t67) * pkin(5) + t36) -m(7) * t36];
taug  = t1(:);
