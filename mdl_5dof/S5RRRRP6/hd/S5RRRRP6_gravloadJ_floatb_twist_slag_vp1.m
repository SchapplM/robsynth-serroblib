% Calculate Gravitation load on the joints for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:49
% EndTime: 2019-12-31 21:52:52
% DurationCPUTime: 0.63s
% Computational Cost: add. (315->109), mult. (392->156), div. (0->0), fcn. (356->8), ass. (0->51)
t78 = rSges(5,3) + pkin(8);
t77 = rSges(6,1) + pkin(4);
t29 = qJ(2) + qJ(3);
t26 = sin(t29);
t27 = cos(t29);
t76 = t27 * rSges(4,1) - t26 * rSges(4,2);
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t75 = g(1) * t36 + g(2) * t33;
t34 = cos(qJ(4));
t24 = t34 * pkin(4) + pkin(3);
t30 = -qJ(5) - pkin(8);
t41 = t27 * t24 + (rSges(6,3) - t30) * t26;
t45 = t27 * pkin(3) + t78 * t26;
t74 = t75 * t26;
t32 = sin(qJ(2));
t73 = pkin(2) * t32;
t31 = sin(qJ(4));
t72 = pkin(4) * t31;
t69 = rSges(3,3) + pkin(6);
t66 = t27 * t31;
t65 = t27 * t33;
t64 = t27 * t34;
t63 = t27 * t36;
t62 = t33 * t31;
t61 = t33 * t34;
t60 = t36 * t31;
t59 = t36 * t34;
t37 = -pkin(7) - pkin(6);
t58 = rSges(4,3) - t37;
t55 = t26 * t62;
t57 = rSges(5,2) * t55 + t78 * t65;
t54 = t26 * t60;
t56 = rSges(5,2) * t54 + t78 * t63;
t53 = -rSges(5,1) * t34 - pkin(3);
t52 = -t37 + t72;
t51 = -rSges(6,1) * t34 - t24;
t35 = cos(qJ(2));
t49 = t35 * rSges(3,1) - t32 * rSges(3,2);
t46 = pkin(1) + t49;
t3 = -t27 * t60 + t61;
t1 = t27 * t62 + t59;
t44 = rSges(5,1) * t64 - rSges(5,2) * t66 + t45;
t42 = g(1) * (rSges(6,2) * t54 + rSges(6,3) * t63) + g(2) * (rSges(6,2) * t55 + rSges(6,3) * t65);
t40 = rSges(6,1) * t64 - rSges(6,2) * t66 + t41;
t28 = t35 * pkin(2);
t25 = t28 + pkin(1);
t10 = t36 * t25;
t4 = t27 * t59 + t62;
t2 = -t27 * t61 + t60;
t5 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - t36 * rSges(2,2)) + g(2) * (t36 * rSges(2,1) - t33 * rSges(2,2))) - m(3) * ((g(1) * t69 + g(2) * t46) * t36 + (-g(1) * t46 + g(2) * t69) * t33) - m(4) * (g(2) * t10 + (g(1) * t58 + g(2) * t76) * t36 + (g(1) * (-t25 - t76) + g(2) * t58) * t33) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t10) + (-g(1) * t37 + g(2) * t45) * t36 + (g(1) * (-t25 - t45) - g(2) * t37) * t33) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10) + (g(1) * t52 + g(2) * t41) * t36 + (g(1) * (-t25 - t41) + g(2) * t52) * t33), -m(3) * (g(3) * t49 + t75 * (-rSges(3,1) * t32 - rSges(3,2) * t35)) - m(4) * (g(3) * (t28 + t76) + t75 * (-rSges(4,1) * t26 - rSges(4,2) * t27 - t73)) - m(5) * (g(1) * (-t36 * t73 + t56) + g(2) * (-t33 * t73 + t57) + g(3) * (t28 + t44) + t53 * t74) - m(6) * (g(3) * (t28 + t40) + t42 + t75 * (t51 * t26 - t27 * t30 - t73)), -m(4) * g(3) * t76 - m(5) * (g(1) * t56 + g(2) * t57 + g(3) * t44) - m(6) * (g(3) * t40 + t42) + t75 * ((m(4) * rSges(4,2) + m(6) * t30) * t27 + (m(4) * rSges(4,1) - m(5) * t53 - m(6) * t51) * t26), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) + t2 * rSges(5,2))) - m(6) * (g(1) * (-t4 * rSges(6,2) + t77 * t3) + g(2) * (t2 * rSges(6,2) - t77 * t1)) + (-m(5) * (-rSges(5,1) * t31 - rSges(5,2) * t34) - m(6) * (-rSges(6,1) * t31 - rSges(6,2) * t34 - t72)) * g(3) * t26, -m(6) * (-g(3) * t27 + t74)];
taug = t5(:);
