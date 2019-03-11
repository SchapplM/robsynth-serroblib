% Calculate Gravitation load on the joints for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:30
% EndTime: 2019-03-09 04:31:31
% DurationCPUTime: 0.66s
% Computational Cost: add. (431->138), mult. (509->176), div. (0->0), fcn. (503->8), ass. (0->48)
t29 = sin(qJ(3));
t59 = g(3) * t29;
t57 = rSges(7,1) + pkin(5);
t27 = qJ(1) + pkin(9);
t23 = cos(t27);
t22 = sin(t27);
t61 = g(2) * t22;
t66 = g(1) * t23 + t61;
t65 = -m(6) - m(7);
t30 = sin(qJ(1));
t64 = pkin(1) * t30;
t63 = g(1) * t22;
t31 = cos(qJ(4));
t60 = t31 * qJ(5) * t59;
t32 = cos(qJ(3));
t25 = t32 * pkin(3);
t58 = -rSges(6,1) - pkin(4);
t56 = t22 * t32;
t55 = t23 * t29;
t54 = t23 * t32;
t28 = sin(qJ(4));
t53 = t28 * t32;
t52 = t31 * t32;
t51 = t29 * pkin(8) + t25;
t50 = rSges(7,2) + qJ(5);
t49 = rSges(6,3) + qJ(5);
t48 = -rSges(7,3) - qJ(6);
t47 = -pkin(4) - t57;
t33 = cos(qJ(1));
t26 = t33 * pkin(1);
t46 = t23 * pkin(2) + t22 * pkin(7) + t26;
t45 = -pkin(2) - t25;
t44 = t23 * pkin(7) - t64;
t43 = t23 * t48;
t42 = pkin(4) * t52 + qJ(5) * t53 + t51;
t41 = pkin(3) * t54 + pkin(8) * t55 + t46;
t40 = rSges(4,1) * t32 - rSges(4,2) * t29;
t9 = t22 * t28 + t23 * t52;
t38 = t9 * pkin(4) + t41;
t6 = t22 * t53 + t23 * t31;
t7 = t22 * t52 - t23 * t28;
t36 = -t7 * pkin(4) - qJ(5) * t6 + t44;
t12 = pkin(8) * t54;
t10 = pkin(8) * t56;
t8 = -t22 * t31 + t23 * t53;
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - rSges(2,2) * t33) + g(2) * (rSges(2,1) * t33 - t30 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t22 - rSges(3,2) * t23 - t64) + g(2) * (rSges(3,1) * t23 - rSges(3,2) * t22 + t26)) - m(4) * (g(1) * (rSges(4,3) * t23 + t44) + g(2) * (rSges(4,1) * t54 - rSges(4,2) * t55 + t46) + (g(1) * (-pkin(2) - t40) + g(2) * rSges(4,3)) * t22) - m(5) * (g(1) * (-rSges(5,1) * t7 + rSges(5,2) * t6 + t44) + g(2) * (rSges(5,1) * t9 - rSges(5,2) * t8 + rSges(5,3) * t55 + t41) + ((-rSges(5,3) - pkin(8)) * t29 + t45) * t63) - m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,3) * t6 + t36) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t55 + t49 * t8 + t38) + ((-rSges(6,2) - pkin(8)) * t29 + t45) * t63) - m(7) * (g(1) * (-rSges(7,2) * t6 - t57 * t7 + t36) + g(2) * (t29 * t43 + t50 * t8 + t57 * t9 + t38) + ((-pkin(8) - t48) * t29 + t45) * t63) (-m(3) - m(4) - m(5) + t65) * g(3), -m(4) * (g(3) * t40 + t66 * (-rSges(4,1) * t29 - rSges(4,2) * t32)) - m(5) * (g(1) * (rSges(5,3) * t54 + t12) + g(2) * (rSges(5,3) * t56 + t10) + g(3) * (rSges(5,1) * t52 - rSges(5,2) * t53 + t51) + (g(3) * rSges(5,3) + t66 * (-rSges(5,1) * t31 + rSges(5,2) * t28 - pkin(3))) * t29) - m(6) * (g(1) * (rSges(6,2) * t54 + t12) + g(2) * (rSges(6,2) * t56 + t10) + g(3) * (rSges(6,1) * t52 + rSges(6,3) * t53 + t42) + (g(3) * rSges(6,2) + t66 * (-t49 * t28 + t58 * t31 - pkin(3))) * t29) - m(7) * (g(1) * t12 + g(2) * t10 + g(3) * t42 + (g(3) * (rSges(7,2) * t28 + t57 * t31) + g(1) * t43 + t48 * t61) * t32 + (g(3) * t48 + t66 * (-t50 * t28 + t47 * t31 - pkin(3))) * t29) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (-rSges(6,1) * t8 + t49 * t9 - t4) + g(2) * (-rSges(6,1) * t6 + t49 * t7 - t2) + t60) - m(7) * (g(1) * (t50 * t9 - t57 * t8 - t4) + g(2) * (t50 * t7 - t57 * t6 - t2) + t60) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t31 + (m(5) * rSges(5,1) - m(6) * t58 - m(7) * t47) * t28) * t59, t65 * (g(1) * t8 + g(2) * t6 + t28 * t59) -m(7) * (g(3) * t32 - t66 * t29)];
taug  = t1(:);
