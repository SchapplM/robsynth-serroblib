% Calculate Gravitation load on the joints for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:40
% EndTime: 2019-12-31 21:59:42
% DurationCPUTime: 0.58s
% Computational Cost: add. (309->120), mult. (440->170), div. (0->0), fcn. (420->8), ass. (0->49)
t36 = -pkin(8) - pkin(7);
t48 = rSges(5,3) - t36;
t47 = rSges(6,3) + qJ(5) - t36;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t64 = g(1) * t35 + g(2) * t32;
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t52 = rSges(4,3) + pkin(7);
t63 = t34 * pkin(2) + t52 * t31;
t29 = qJ(3) + qJ(4);
t22 = sin(t29);
t23 = cos(t29);
t50 = t32 * t34;
t10 = t22 * t35 - t23 * t50;
t9 = t22 * t50 + t23 * t35;
t62 = -t9 * rSges(5,1) + t10 * rSges(5,2);
t61 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t49 = t34 * t35;
t11 = -t22 * t49 + t23 * t32;
t12 = t22 * t32 + t23 * t49;
t60 = t11 * rSges(6,1) - t12 * rSges(6,2);
t59 = t11 * rSges(5,1) - t12 * rSges(5,2);
t30 = sin(qJ(3));
t58 = pkin(3) * t30;
t57 = pkin(4) * t22;
t54 = g(3) * t31;
t51 = rSges(3,2) * t31;
t33 = cos(qJ(3));
t25 = t33 * pkin(3);
t19 = pkin(4) * t23 + t25;
t46 = t35 * pkin(1) + t32 * pkin(6);
t45 = rSges(3,1) * t34 - t51;
t43 = -rSges(5,1) * t22 - rSges(5,2) * t23;
t42 = -rSges(6,1) * t22 - rSges(6,2) * t23;
t41 = rSges(4,1) * t33 - rSges(4,2) * t30 + pkin(2);
t15 = -t30 * t49 + t32 * t33;
t13 = t30 * t50 + t33 * t35;
t21 = t25 + pkin(2);
t40 = rSges(5,1) * t23 - rSges(5,2) * t22 + t21;
t17 = pkin(2) + t19;
t39 = rSges(6,1) * t23 - rSges(6,2) * t22 + t17;
t38 = t34 * t21 + t31 * t48;
t37 = t17 * t34 + t31 * t47;
t26 = t35 * pkin(6);
t18 = t57 + t58;
t16 = t30 * t32 + t33 * t49;
t14 = t30 * t35 - t33 * t50;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t32 - rSges(2,2) * t35) + g(2) * (rSges(2,1) * t35 - rSges(2,2) * t32)) - m(3) * (g(1) * (rSges(3,3) * t35 + t26) + g(2) * (rSges(3,1) * t49 - t35 * t51 + t46) + (g(1) * (-pkin(1) - t45) + g(2) * rSges(3,3)) * t32) - m(4) * ((rSges(4,1) * t16 + rSges(4,2) * t15 + t35 * t63 + t46) * g(2) + (rSges(4,1) * t14 + rSges(4,2) * t13 + t26 + (-pkin(1) - t63) * t32) * g(1)) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t26) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t46) + (g(1) * t58 + g(2) * t38) * t35 + (g(1) * (-pkin(1) - t38) + g(2) * t58) * t32) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 + t26) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t11 + t46) + (g(1) * t18 + g(2) * t37) * t35 + (g(1) * (-pkin(1) - t37) + g(2) * t18) * t32), -m(3) * (g(3) * t45 + t64 * (-rSges(3,1) * t31 - rSges(3,2) * t34)) - m(4) * ((g(3) * t41 + t52 * t64) * t34 + (g(3) * t52 - t41 * t64) * t31) - m(5) * ((g(3) * t40 + t48 * t64) * t34 + (g(3) * t48 - t40 * t64) * t31) - m(6) * ((g(3) * t39 + t47 * t64) * t34 + (g(3) * t47 - t39 * t64) * t31), -m(4) * (g(1) * (rSges(4,1) * t15 - rSges(4,2) * t16) + g(2) * (-rSges(4,1) * t13 + rSges(4,2) * t14)) - m(5) * (g(1) * (t15 * pkin(3) + t59) + g(2) * (-t13 * pkin(3) + t62)) - m(6) * (g(1) * (-t18 * t49 + t19 * t32 + t60) + g(2) * (-t18 * t50 - t19 * t35 + t61)) + (-m(4) * (-rSges(4,1) * t30 - rSges(4,2) * t33) - m(5) * (t43 - t58) - m(6) * (-t18 + t42)) * t54, -m(5) * (g(1) * t59 + g(2) * t62) - m(6) * (g(1) * (t11 * pkin(4) + t60) + g(2) * (-t9 * pkin(4) + t61)) + (-m(5) * t43 - m(6) * (t42 - t57)) * t54, -m(6) * (-g(3) * t34 + t31 * t64)];
taug = t1(:);
