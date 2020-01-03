% Calculate Gravitation load on the joints for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:52
% EndTime: 2019-12-31 20:54:54
% DurationCPUTime: 0.42s
% Computational Cost: add. (280->79), mult. (254->96), div. (0->0), fcn. (208->8), ass. (0->45)
t67 = rSges(6,3) + qJ(5);
t66 = rSges(6,1) + pkin(4);
t26 = qJ(2) + qJ(3);
t21 = pkin(8) + t26;
t18 = cos(t21);
t65 = t67 * t18;
t30 = cos(qJ(1));
t64 = g(2) * t30;
t17 = sin(t21);
t38 = t18 * rSges(5,1) - t17 * rSges(5,2);
t22 = sin(t26);
t23 = cos(t26);
t44 = t23 * rSges(4,1) - t22 * rSges(4,2);
t28 = sin(qJ(1));
t62 = g(1) * t30 + g(2) * t28;
t61 = t62 * t17;
t33 = t67 * t17 + t66 * t18;
t19 = pkin(3) * t23;
t29 = cos(qJ(2));
t24 = t29 * pkin(2);
t20 = t24 + pkin(1);
t4 = t19 + t20;
t60 = t4 * t64;
t31 = -pkin(7) - pkin(6);
t59 = t65 * t28;
t58 = t65 * t30;
t57 = pkin(3) * t22;
t27 = sin(qJ(2));
t54 = t27 * pkin(2);
t52 = rSges(3,3) + pkin(6);
t25 = -qJ(4) + t31;
t48 = rSges(6,2) - t25;
t47 = rSges(4,3) - t31;
t46 = rSges(5,3) - t25;
t43 = t19 + t38;
t42 = t19 + t33;
t41 = t29 * rSges(3,1) - t27 * rSges(3,2);
t39 = -rSges(4,1) * t22 - rSges(4,2) * t23;
t37 = -rSges(5,1) * t17 - rSges(5,2) * t18;
t36 = pkin(1) + t41;
t35 = t20 + t44;
t5 = -t54 - t57;
t3 = t30 * t5;
t2 = t28 * t5;
t1 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - t30 * rSges(2,2)) + g(2) * (t30 * rSges(2,1) - t28 * rSges(2,2))) - m(3) * ((g(1) * t52 + g(2) * t36) * t30 + (-g(1) * t36 + g(2) * t52) * t28) - m(4) * ((g(1) * t47 + g(2) * t35) * t30 + (-g(1) * t35 + g(2) * t47) * t28) - m(5) * (t60 + (g(1) * t46 + g(2) * t38) * t30 + (g(1) * (-t38 - t4) + g(2) * t46) * t28) - m(6) * (t60 + (g(1) * t48 + g(2) * t33) * t30 + (g(1) * (-t33 - t4) + g(2) * t48) * t28), -m(3) * (g(3) * t41 + t62 * (-rSges(3,1) * t27 - rSges(3,2) * t29)) - m(4) * (g(3) * (t24 + t44) + t62 * (t39 - t54)) - m(5) * (g(1) * (t37 * t30 + t3) + g(2) * (t37 * t28 + t2) + g(3) * (t24 + t43)) - m(6) * (g(1) * (t3 + t58) + g(2) * (t2 + t59) + g(3) * (t24 + t42) - t66 * t61), -m(6) * (g(1) * t58 + g(2) * t59) + (-m(4) * t44 - m(5) * t43 - m(6) * t42) * g(3) + t62 * (-m(4) * t39 - m(5) * (t37 - t57) - m(6) * (-t17 * t66 - t57)), (-m(5) - m(6)) * (g(1) * t28 - t64), -m(6) * (-g(3) * t18 + t61)];
taug = t1(:);
