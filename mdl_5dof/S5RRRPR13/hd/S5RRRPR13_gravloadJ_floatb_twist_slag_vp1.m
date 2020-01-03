% Calculate Gravitation load on the joints for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:18
% EndTime: 2019-12-31 21:43:21
% DurationCPUTime: 0.91s
% Computational Cost: add. (405->146), mult. (969->216), div. (0->0), fcn. (1144->10), ass. (0->62)
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t90 = -rSges(6,1) * t42 - rSges(6,2) * t46;
t47 = cos(qJ(3));
t43 = sin(qJ(3));
t65 = qJ(4) * t43;
t87 = -pkin(3) * t47 - t65;
t83 = pkin(4) + pkin(8);
t84 = rSges(6,1) * t46 - rSges(6,2) * t42 + t83;
t78 = pkin(9) + rSges(6,3);
t85 = t90 * t43;
t41 = sin(pkin(5));
t81 = g(3) * t41;
t80 = rSges(5,1) + pkin(8);
t79 = rSges(4,3) + pkin(8);
t77 = cos(qJ(1));
t44 = sin(qJ(2));
t76 = t41 * t44;
t45 = sin(qJ(1));
t75 = t41 * t45;
t74 = t41 * t47;
t48 = cos(qJ(2));
t73 = t41 * t48;
t71 = t42 * t48;
t69 = t46 * t48;
t68 = t47 * t48;
t67 = t77 * pkin(1) + pkin(7) * t75;
t66 = pkin(2) * t73 + pkin(8) * t76;
t64 = rSges(5,3) + qJ(4);
t63 = cos(pkin(5));
t54 = t63 * t77;
t25 = t44 * t45 - t48 * t54;
t19 = t25 * pkin(2);
t62 = t87 * t25 - t19;
t57 = t45 * t63;
t27 = t77 * t44 + t48 * t57;
t21 = t27 * pkin(2);
t61 = t87 * t27 - t21;
t28 = -t44 * t57 + t77 * t48;
t60 = t28 * pkin(2) + t67;
t59 = t41 * t77;
t58 = -t45 * pkin(1) + pkin(7) * t59;
t26 = t44 * t54 + t45 * t48;
t9 = t26 * t47 - t43 * t59;
t13 = t28 * t47 + t43 * t75;
t56 = t13 * pkin(3) + t60;
t55 = -t26 * pkin(2) + t58;
t53 = rSges(4,1) * t47 - rSges(4,2) * t43;
t52 = rSges(5,2) * t47 - rSges(5,3) * t43;
t51 = -pkin(3) * t9 + t55;
t50 = g(3) * (t41 * pkin(3) * t68 + t65 * t73 + t66);
t49 = qJ(4) - t90;
t8 = t26 * t43 + t47 * t59;
t24 = t63 * t43 + t44 * t74;
t23 = t43 * t76 - t63 * t47;
t18 = t23 * pkin(3);
t12 = t28 * t43 - t45 * t74;
t6 = t12 * pkin(3);
t4 = t8 * pkin(3);
t3 = t12 * t42 + t27 * t46;
t2 = t12 * t46 - t27 * t42;
t1 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) - t45 * rSges(2,2))) - m(3) * (g(1) * (-t26 * rSges(3,1) + t25 * rSges(3,2) + rSges(3,3) * t59 + t58) + g(2) * (rSges(3,1) * t28 - rSges(3,2) * t27 + rSges(3,3) * t75 + t67)) - m(4) * (g(1) * (-rSges(4,1) * t9 + rSges(4,2) * t8 - t79 * t25 + t55) + g(2) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t79 * t27 + t60)) - m(5) * (g(1) * (rSges(5,2) * t9 - t80 * t25 - t64 * t8 + t51) + g(2) * (-t13 * rSges(5,2) + t64 * t12 + t80 * t27 + t56)) - m(6) * (g(1) * (-t84 * t25 - t49 * t8 - t78 * t9 + t51) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t12 * qJ(4) + t78 * t13 + t83 * t27 + t56)), -m(3) * (g(1) * (-rSges(3,1) * t27 - rSges(3,2) * t28) + g(2) * (-rSges(3,1) * t25 - rSges(3,2) * t26) + (rSges(3,1) * t48 - rSges(3,2) * t44) * t81) - m(4) * (g(1) * (-t53 * t27 + t79 * t28 - t21) + g(2) * (-t53 * t25 + t79 * t26 - t19) + g(3) * t66 + (rSges(4,3) * t44 + t53 * t48) * t81) - m(5) * (g(1) * (t52 * t27 + t80 * t28 + t61) + g(2) * (t52 * t25 + t80 * t26 + t62) + t50 + (rSges(5,1) * t44 - t52 * t48) * t81) - m(6) * (g(1) * (t85 * t27 + t84 * t28 + t61) + g(2) * (t85 * t25 + t84 * t26 + t62) + t50 + (t44 * pkin(4) + (t43 * t71 + t44 * t46) * rSges(6,1) + (-t42 * t44 + t43 * t69) * rSges(6,2)) * t81 + ((-g(1) * t27 - g(2) * t25) * t47 + t68 * t81) * t78), -m(4) * (g(1) * (-rSges(4,1) * t12 - rSges(4,2) * t13) + g(2) * (-rSges(4,1) * t8 - rSges(4,2) * t9) + g(3) * (-rSges(4,1) * t23 - rSges(4,2) * t24)) - m(5) * (g(1) * (t12 * rSges(5,2) + t64 * t13 - t6) + g(2) * (t8 * rSges(5,2) + t64 * t9 - t4) + g(3) * (t23 * rSges(5,2) + t64 * t24 - t18)) + (-g(1) * (-t78 * t12 - t6) - g(2) * (-t78 * t8 - t4) - g(3) * (-t78 * t23 - t18) - (g(1) * t13 + g(2) * t9 + g(3) * t24) * t49) * m(6), (-m(5) - m(6)) * (g(1) * t12 + g(2) * t8 + g(3) * t23), -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * ((-t25 * t42 + t46 * t8) * rSges(6,1) + (-t25 * t46 - t42 * t8) * rSges(6,2)) + g(3) * ((t23 * t46 + t41 * t71) * rSges(6,1) + (-t23 * t42 + t41 * t69) * rSges(6,2)))];
taug = t1(:);
