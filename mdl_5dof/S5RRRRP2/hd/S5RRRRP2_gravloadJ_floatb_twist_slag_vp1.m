% Calculate Gravitation load on the joints for
% S5RRRRP2
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:33
% EndTime: 2019-12-05 18:47:34
% DurationCPUTime: 0.33s
% Computational Cost: add. (290->84), mult. (234->108), div. (0->0), fcn. (186->8), ass. (0->50)
t70 = rSges(4,3) + pkin(7);
t36 = -pkin(8) - pkin(7);
t69 = rSges(5,3) - t36;
t68 = rSges(6,3) + qJ(5) - t36;
t34 = cos(qJ(3));
t58 = rSges(4,1) * t34;
t67 = -pkin(2) - t58;
t30 = qJ(3) + qJ(4);
t26 = cos(t30);
t66 = (-rSges(6,1) - pkin(4)) * t26;
t31 = qJ(1) + qJ(2);
t25 = sin(t31);
t54 = t25 * t26;
t24 = sin(t30);
t56 = t24 * t25;
t65 = rSges(6,1) * t56 + rSges(6,2) * t54;
t64 = rSges(5,1) * t56 + rSges(5,2) * t54;
t32 = sin(qJ(3));
t63 = pkin(3) * t32;
t62 = pkin(4) * t24;
t27 = cos(t31);
t61 = g(3) * t27;
t33 = sin(qJ(1));
t60 = t33 * pkin(1);
t35 = cos(qJ(1));
t59 = t35 * pkin(1);
t28 = t34 * pkin(3);
t23 = t28 + pkin(2);
t57 = rSges(4,2) * t32;
t55 = t24 * t27;
t53 = t25 * t32;
t17 = t26 * rSges(5,1);
t52 = rSges(4,2) * t53 + t70 * t27;
t51 = -t23 + t66;
t50 = -rSges(3,1) * t27 + t25 * rSges(3,2);
t49 = -t23 - t17;
t48 = -rSges(5,2) * t24 + t17;
t47 = -rSges(6,2) * t24 - t66;
t46 = -rSges(3,1) * t25 - rSges(3,2) * t27;
t45 = rSges(4,1) * t32 + rSges(4,2) * t34;
t44 = -rSges(5,1) * t24 - rSges(5,2) * t26;
t43 = -rSges(6,1) * t24 - rSges(6,2) * t26;
t42 = (t57 + t67) * t27;
t41 = rSges(6,2) * t55 - t68 * t25 + t51 * t27;
t40 = rSges(6,2) * t56 + t51 * t25 + t68 * t27;
t39 = rSges(5,2) * t56 + t49 * t25 + t69 * t27;
t38 = rSges(5,2) * t55 - t69 * t25 + t49 * t27;
t37 = (-g(2) * t70 + g(3) * t67) * t25;
t2 = -t62 - t63;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t35 + rSges(2,2) * t33) + g(3) * (-rSges(2,1) * t33 - rSges(2,2) * t35)) - m(3) * (g(2) * (t50 - t59) + g(3) * (t46 - t60)) - m(4) * (g(2) * (t42 - t59) + g(3) * (t52 - t60) + t37) - m(5) * (g(2) * (t38 - t59) + g(3) * (t39 - t60)) - m(6) * (g(2) * (t41 - t59) + g(3) * (t40 - t60)), -m(3) * (g(2) * t50 + g(3) * t46) - m(4) * (g(2) * t42 + g(3) * t52 + t37) - m(5) * (g(2) * t38 + g(3) * t39) - m(6) * (g(2) * t41 + g(3) * t40), -m(4) * (g(1) * (-t57 + t58) + g(2) * t45 * t25) - m(5) * (g(1) * (t28 + t48) + g(2) * (pkin(3) * t53 + t64)) - m(6) * (g(1) * (t28 + t47) + g(2) * (-t2 * t25 + t65)) + (m(4) * t45 - m(5) * (t44 - t63) - m(6) * (t2 + t43)) * t61, -m(5) * (g(1) * t48 + g(2) * t64) - m(6) * (g(1) * t47 + g(2) * (pkin(4) * t56 + t65)) + (-m(5) * t44 - m(6) * (t43 - t62)) * t61, -m(6) * (g(2) * t27 + g(3) * t25)];
taug = t1(:);
