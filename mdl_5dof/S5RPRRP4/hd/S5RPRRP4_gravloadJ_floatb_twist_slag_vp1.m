% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:30
% EndTime: 2019-12-05 18:05:32
% DurationCPUTime: 0.46s
% Computational Cost: add. (233->103), mult. (324->139), div. (0->0), fcn. (316->8), ass. (0->46)
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t61 = -pkin(2) * t32 - pkin(1) + (-rSges(4,3) - pkin(6)) * t31;
t37 = -pkin(7) - pkin(6);
t30 = qJ(3) + qJ(4);
t25 = sin(t30);
t36 = cos(qJ(1));
t47 = t36 * t25;
t26 = cos(t30);
t34 = sin(qJ(1));
t50 = t34 * t26;
t10 = t32 * t50 - t47;
t46 = t36 * t26;
t51 = t34 * t25;
t9 = t32 * t51 + t46;
t60 = t9 * rSges(6,1) + t10 * rSges(6,2);
t59 = t9 * rSges(5,1) + t10 * rSges(5,2);
t11 = t32 * t47 - t50;
t12 = -t32 * t46 - t51;
t58 = -t11 * rSges(5,1) + t12 * rSges(5,2);
t57 = -t11 * rSges(6,1) + t12 * rSges(6,2);
t56 = pkin(4) * t25;
t55 = g(1) * t31;
t54 = g(2) * t36;
t33 = sin(qJ(3));
t53 = t33 * pkin(3);
t18 = t53 + t56;
t52 = t18 * t32;
t49 = t34 * t33;
t35 = cos(qJ(3));
t48 = t34 * t35;
t45 = t36 * t33;
t44 = t36 * t35;
t28 = t35 * pkin(3);
t19 = pkin(4) * t26 + t28;
t42 = -rSges(5,1) * t25 - rSges(5,2) * t26;
t41 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t40 = -rSges(3,1) * t32 + rSges(3,2) * t31 - pkin(1);
t15 = t32 * t45 - t48;
t13 = t32 * t49 + t44;
t39 = -(t28 + pkin(2)) * t32 - pkin(1) + (-rSges(5,3) + t37) * t31;
t38 = -(pkin(2) + t19) * t32 - pkin(1) + (-rSges(6,3) - qJ(5) + t37) * t31;
t27 = t36 * qJ(2);
t16 = -t32 * t44 - t49;
t14 = t32 * t48 - t45;
t1 = [-m(2) * (g(2) * (-t36 * rSges(2,1) + t34 * rSges(2,2)) + g(3) * (-t34 * rSges(2,1) - t36 * rSges(2,2))) - m(3) * (g(3) * t27 + (g(3) * rSges(3,3) + g(2) * t40) * t36 + (g(2) * (-rSges(3,3) - qJ(2)) + g(3) * t40) * t34) - m(4) * (g(2) * (t16 * rSges(4,1) + t15 * rSges(4,2)) + g(3) * (-t14 * rSges(4,1) + t13 * rSges(4,2) + t27) + t61 * t54 + (-g(2) * qJ(2) + g(3) * t61) * t34) - m(5) * (g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2)) + g(3) * (-t10 * rSges(5,1) + t9 * rSges(5,2) + t27) + (g(2) * t39 + g(3) * t53) * t36 + (g(2) * (-qJ(2) - t53) + g(3) * t39) * t34) - m(6) * (g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2)) + g(3) * (-t10 * rSges(6,1) + t9 * rSges(6,2) + t27) + (g(2) * t38 + g(3) * t18) * t36 + (g(2) * (-qJ(2) - t18) + g(3) * t38) * t34), (-m(3) - m(4) - m(5) - m(6)) * (g(3) * t34 + t54), -m(4) * (g(2) * (t13 * rSges(4,1) + t14 * rSges(4,2)) + g(3) * (-t15 * rSges(4,1) + t16 * rSges(4,2))) - m(5) * (g(2) * (t13 * pkin(3) + t59) + g(3) * (-t15 * pkin(3) + t58)) - m(6) * (g(2) * (t36 * t19 + t34 * t52 + t60) + g(3) * (t34 * t19 - t36 * t52 + t57)) + (-m(4) * (-rSges(4,1) * t33 - rSges(4,2) * t35) - m(5) * (t42 - t53) - m(6) * (-t18 + t41)) * t55, -m(5) * (g(2) * t59 + g(3) * t58) - m(6) * (g(2) * (t9 * pkin(4) + t60) + g(3) * (-t11 * pkin(4) + t57)) + (-m(5) * t42 - m(6) * (t41 - t56)) * t55, -m(6) * (-g(1) * t32 + (-g(2) * t34 + g(3) * t36) * t31)];
taug = t1(:);
