% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:04
% DurationCPUTime: 0.25s
% Computational Cost: add. (349->67), mult. (217->87), div. (0->0), fcn. (169->10), ass. (0->44)
t39 = cos(qJ(4));
t35 = qJ(4) + qJ(5);
t27 = sin(t35);
t29 = cos(t35);
t53 = t29 * rSges(6,1) - t27 * rSges(6,2);
t73 = t39 * pkin(4) + t53;
t72 = rSges(5,3) + pkin(8);
t71 = rSges(6,3) + pkin(9) + pkin(8);
t37 = sin(qJ(4));
t70 = t39 * rSges(5,1) - t37 * rSges(5,2);
t69 = rSges(6,1) * t27 + rSges(6,2) * t29;
t68 = pkin(3) + t70;
t67 = pkin(3) + t73;
t36 = qJ(1) + qJ(2);
t31 = qJ(3) + t36;
t25 = cos(t31);
t66 = t69 * t25;
t65 = pkin(4) * t37;
t24 = sin(t31);
t64 = g(2) * t24;
t58 = t24 * rSges(4,1) + t25 * rSges(4,2);
t28 = sin(t36);
t30 = cos(t36);
t57 = t28 * rSges(3,1) + t30 * rSges(3,2);
t22 = pkin(2) * t28;
t56 = t22 + t58;
t55 = t30 * rSges(3,1) - rSges(3,2) * t28;
t54 = t25 * rSges(4,1) - rSges(4,2) * t24;
t23 = pkin(2) * t30;
t52 = t23 + t54;
t51 = rSges(5,1) * t37 + rSges(5,2) * t39;
t49 = t72 * t24 + t68 * t25;
t48 = t23 + t49;
t47 = t67 * t24 - t71 * t25;
t46 = t71 * t24 + t67 * t25;
t45 = t68 * t24 - t72 * t25;
t44 = t22 + t47;
t43 = t23 + t46;
t42 = t22 + t45;
t40 = cos(qJ(1));
t38 = sin(qJ(1));
t34 = t40 * pkin(1);
t32 = t38 * pkin(1);
t1 = [-m(2) * (g(2) * (t40 * rSges(2,1) - t38 * rSges(2,2)) + g(3) * (t38 * rSges(2,1) + t40 * rSges(2,2))) - m(3) * (g(2) * (t34 + t55) + g(3) * (t32 + t57)) - m(4) * (g(2) * (t34 + t52) + g(3) * (t32 + t56)) - m(5) * (g(2) * (t34 + t48) + g(3) * (t32 + t42)) - m(6) * (g(2) * (t34 + t43) + g(3) * (t32 + t44)), -m(3) * (g(2) * t55 + g(3) * t57) - m(4) * (g(2) * t52 + g(3) * t56) - m(5) * (g(2) * t48 + g(3) * t42) - m(6) * (g(2) * t43 + g(3) * t44), -m(4) * (g(2) * t54 + g(3) * t58) - m(5) * (g(2) * t49 + g(3) * t45) - m(6) * (g(2) * t46 + g(3) * t47), -m(5) * (g(3) * t51 * t25 + g(1) * t70) - m(6) * (g(1) * t73 + g(3) * (t25 * t65 + t66)) + (m(5) * t51 - m(6) * (-t69 - t65)) * t64, -m(6) * (g(1) * t53 + g(3) * t66 - t64 * t69)];
taug = t1(:);
