% Calculate Gravitation load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:53
% EndTime: 2019-12-05 18:23:55
% DurationCPUTime: 0.46s
% Computational Cost: add. (183->76), mult. (255->112), div. (0->0), fcn. (217->8), ass. (0->43)
t20 = qJ(2) + qJ(4);
t18 = sin(t20);
t19 = cos(t20);
t22 = sin(qJ(5));
t54 = rSges(6,2) * t22;
t68 = rSges(6,3) + pkin(4);
t69 = t18 * t54 + t19 * t68;
t26 = cos(qJ(2));
t28 = pkin(2) + pkin(1);
t17 = t26 * t28;
t64 = t19 * rSges(5,1) - t18 * rSges(5,2);
t67 = -t64 - t17;
t65 = t68 * t18;
t23 = sin(qJ(2));
t57 = rSges(4,1) + pkin(1);
t30 = -t23 * rSges(4,2) + t57 * t26;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t63 = -g(1) * t24 + g(2) * t27;
t62 = g(1) * t27 + g(2) * t24;
t25 = cos(qJ(5));
t55 = rSges(6,1) * t25;
t49 = t23 * t28;
t48 = t24 * t22;
t47 = t24 * t25;
t46 = t27 * t22;
t45 = t27 * t25;
t21 = -pkin(3) - qJ(3);
t44 = rSges(5,3) - t21;
t43 = rSges(4,3) + qJ(3);
t42 = t69 * t24;
t41 = t69 * t27;
t40 = t18 * t55;
t37 = t26 * rSges(3,1) - t23 * rSges(3,2);
t34 = -rSges(5,1) * t18 - rSges(5,2) * t19;
t33 = -t40 - t49;
t32 = (-t54 + t55) * t19 + t65;
t13 = t27 * t17;
t4 = t19 * t45 + t48;
t3 = -t19 * t46 + t47;
t2 = -t19 * t47 + t46;
t1 = t19 * t48 + t45;
t5 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(1) * (t27 * rSges(3,3) - t37 * t24) + g(2) * (t24 * rSges(3,3) + t37 * t27)) - m(4) * ((g(1) * t43 + g(2) * t30) * t27 + (-g(1) * t30 + g(2) * t43) * t24) - m(5) * (g(2) * t13 + (g(1) * t44 + g(2) * t64) * t27 + (g(1) * t67 + g(2) * t44) * t24) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t24 * t17 - t27 * t21) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) - t24 * t21 + t13) + t63 * t65), -m(3) * (g(3) * t37 + t62 * (-rSges(3,1) * t23 - rSges(3,2) * t26)) - m(4) * (g(3) * t30 + t62 * (-rSges(4,2) * t26 - t57 * t23)) - m(5) * (-g(3) * t67 + t62 * (t34 - t49)) - m(6) * (g(1) * (t33 * t27 + t41) + g(2) * (t33 * t24 + t42) + g(3) * (t17 + t32)), -(-m(4) - m(5) - m(6)) * t63, -m(5) * (g(3) * t64 + t62 * t34) - m(6) * (g(1) * (-t27 * t40 + t41) + g(2) * (-t24 * t40 + t42) + g(3) * t32), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t22 - rSges(6,2) * t25) * t18)];
taug = t5(:);
