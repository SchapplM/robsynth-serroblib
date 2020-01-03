% Calculate Gravitation load on the joints for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:36
% EndTime: 2019-12-31 18:48:37
% DurationCPUTime: 0.35s
% Computational Cost: add. (245->64), mult. (211->76), div. (0->0), fcn. (171->8), ass. (0->33)
t60 = rSges(6,3) + qJ(5);
t59 = rSges(6,1) + pkin(4);
t26 = sin(qJ(1));
t27 = cos(qJ(1));
t52 = g(1) * t27 + g(2) * t26;
t55 = g(2) * t27;
t22 = pkin(8) + qJ(3);
t19 = qJ(4) + t22;
t14 = sin(t19);
t15 = cos(t19);
t53 = t15 * rSges(5,1) - rSges(5,2) * t14;
t51 = t52 * t14;
t50 = t60 * t14 + t59 * t15;
t18 = cos(t22);
t13 = pkin(3) * t18;
t24 = cos(pkin(8));
t16 = t24 * pkin(2) + pkin(1);
t2 = t13 + t16;
t49 = t2 * t55;
t17 = sin(t22);
t48 = pkin(3) * t17;
t25 = -pkin(6) - qJ(2);
t21 = -pkin(7) + t25;
t42 = rSges(6,2) - t21;
t41 = rSges(4,3) - t25;
t40 = rSges(5,3) - t21;
t38 = rSges(3,3) + qJ(2);
t35 = rSges(4,1) * t18 - t17 * rSges(4,2);
t33 = -rSges(5,1) * t14 - rSges(5,2) * t15;
t32 = rSges(3,1) * t24 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t31 = t16 + t35;
t30 = t52 * t60 * t15;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t26 * rSges(2,2))) - m(3) * ((g(1) * t38 + g(2) * t32) * t27 + (-g(1) * t32 + g(2) * t38) * t26) - m(4) * ((g(1) * t41 + g(2) * t31) * t27 + (-g(1) * t31 + g(2) * t41) * t26) - m(5) * (t49 + (g(1) * t40 + g(2) * t53) * t27 + (g(1) * (-t2 - t53) + g(2) * t40) * t26) - m(6) * (t49 + (g(1) * t42 + g(2) * t50) * t27 + (g(1) * (-t2 - t50) + g(2) * t42) * t26), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t26 - t55), -m(6) * t30 + (-m(4) * t35 - m(5) * (t13 + t53) - m(6) * (t13 + t50)) * g(3) + t52 * (-m(4) * (-rSges(4,1) * t17 - rSges(4,2) * t18) - m(5) * (t33 - t48) - m(6) * (-t14 * t59 - t48)), -m(5) * (g(3) * t53 + t52 * t33) - m(6) * (g(3) * t50 - t51 * t59 + t30), -m(6) * (-g(3) * t15 + t51)];
taug = t1(:);
