% Calculate Gravitation load on the joints for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:30
% EndTime: 2019-12-31 21:07:33
% DurationCPUTime: 0.64s
% Computational Cost: add. (225->127), mult. (492->167), div. (0->0), fcn. (496->6), ass. (0->44)
t25 = sin(qJ(2));
t56 = g(3) * t25;
t27 = cos(qJ(3));
t62 = t27 * t56;
t29 = cos(qJ(1));
t26 = sin(qJ(1));
t58 = g(2) * t26;
t61 = g(1) * t29 + t58;
t42 = rSges(6,3) + qJ(5);
t60 = g(1) * t26;
t57 = qJ(4) * t62;
t28 = cos(qJ(2));
t21 = t28 * pkin(2);
t55 = rSges(6,1) + pkin(4);
t54 = rSges(5,2) - pkin(3);
t24 = sin(qJ(3));
t53 = t24 * t28;
t51 = t25 * t29;
t50 = t26 * t28;
t49 = t27 * t28;
t48 = t28 * t29;
t47 = t29 * t24;
t46 = t25 * pkin(7) + t21;
t45 = t29 * pkin(1) + t26 * pkin(6);
t44 = rSges(6,2) + qJ(4);
t43 = rSges(5,3) + qJ(4);
t41 = -pkin(3) - t42;
t40 = -pkin(1) - t21;
t39 = t55 * t29;
t38 = pkin(3) * t49 + qJ(4) * t53 + t46;
t37 = pkin(2) * t48 + pkin(7) * t51 + t45;
t22 = t29 * pkin(6);
t6 = t24 * t50 + t27 * t29;
t7 = t26 * t49 - t47;
t36 = -t7 * pkin(3) - qJ(4) * t6 + t22;
t9 = t26 * t24 + t27 * t48;
t35 = t9 * pkin(3) + t37;
t34 = rSges(3,1) * t28 - rSges(3,2) * t25;
t15 = pkin(7) * t48;
t12 = pkin(7) * t50;
t8 = -t26 * t27 + t28 * t47;
t4 = t8 * pkin(3);
t2 = t6 * pkin(3);
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - rSges(2,2) * t29) + g(2) * (rSges(2,1) * t29 - t26 * rSges(2,2))) - m(3) * (g(1) * (t29 * rSges(3,3) + t22) + g(2) * (rSges(3,1) * t48 - rSges(3,2) * t51 + t45) + (g(1) * (-pkin(1) - t34) + g(2) * rSges(3,3)) * t26) - m(4) * (g(1) * (-rSges(4,1) * t7 + rSges(4,2) * t6 + t22) + g(2) * (t9 * rSges(4,1) - t8 * rSges(4,2) + rSges(4,3) * t51 + t37) + ((-rSges(4,3) - pkin(7)) * t25 + t40) * t60) - m(5) * (g(1) * (rSges(5,2) * t7 - rSges(5,3) * t6 + t36) + g(2) * (rSges(5,1) * t51 - t9 * rSges(5,2) + t43 * t8 + t35) + ((-rSges(5,1) - pkin(7)) * t25 + t40) * t60) - m(6) * (g(1) * (-rSges(6,2) * t6 - t42 * t7 + t36) + g(2) * (t25 * t39 + t42 * t9 + t44 * t8 + t35) + ((-pkin(7) - t55) * t25 + t40) * t60), -m(3) * (g(3) * t34 + t61 * (-rSges(3,1) * t25 - rSges(3,2) * t28)) - m(4) * (g(1) * (rSges(4,3) * t48 + t15) + g(2) * (rSges(4,3) * t50 + t12) + g(3) * (rSges(4,1) * t49 - rSges(4,2) * t53 + t46) + (g(3) * rSges(4,3) + t61 * (-rSges(4,1) * t27 + rSges(4,2) * t24 - pkin(2))) * t25) - m(5) * (g(1) * (rSges(5,1) * t48 + t15) + g(2) * (rSges(5,1) * t50 + t12) + g(3) * (-rSges(5,2) * t49 + rSges(5,3) * t53 + t38) + (g(3) * rSges(5,1) + t61 * (-t43 * t24 + t54 * t27 - pkin(2))) * t25) - m(6) * (g(1) * t15 + g(2) * t12 + g(3) * t38 + (g(3) * (rSges(6,2) * t24 + t42 * t27) + g(1) * t39 + t55 * t58) * t28 + (g(3) * t55 + t61 * (-t44 * t24 + t41 * t27 - pkin(2))) * t25), -m(4) * (g(1) * (-rSges(4,1) * t8 - rSges(4,2) * t9) + g(2) * (-rSges(4,1) * t6 - rSges(4,2) * t7)) - m(5) * (g(1) * (t8 * rSges(5,2) + t43 * t9 - t4) + g(2) * (t6 * rSges(5,2) + t43 * t7 - t2) + t57) - m(6) * (g(1) * (-t42 * t8 + t44 * t9 - t4) + g(2) * (-t42 * t6 + t44 * t7 - t2) + t57) + ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t27 + (m(4) * rSges(4,1) - m(5) * t54 - m(6) * t41) * t24) * t56, (-m(5) - m(6)) * (g(1) * t8 + g(2) * t6 + t24 * t56), -m(6) * (g(1) * t9 + g(2) * t7 + t62)];
taug = t1(:);
