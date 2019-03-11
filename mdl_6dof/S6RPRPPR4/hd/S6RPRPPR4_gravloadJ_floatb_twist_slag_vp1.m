% Calculate Gravitation load on the joints for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:47
% EndTime: 2019-03-09 02:46:48
% DurationCPUTime: 0.70s
% Computational Cost: add. (395->133), mult. (509->176), div. (0->0), fcn. (515->10), ass. (0->56)
t67 = -pkin(8) - rSges(7,3);
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t73 = g(1) * t35 + g(2) * t33;
t72 = -m(6) - m(7);
t31 = -pkin(7) - qJ(2);
t70 = g(2) * t31;
t26 = pkin(9) + qJ(3);
t24 = sin(t26);
t68 = g(3) * t24;
t25 = cos(t26);
t20 = t25 * pkin(3);
t66 = t24 * t35;
t27 = sin(pkin(10));
t65 = t25 * t27;
t29 = cos(pkin(10));
t64 = t25 * t29;
t63 = t25 * t33;
t62 = t25 * t35;
t61 = t29 * t35;
t60 = t33 * t27;
t59 = t33 * t29;
t58 = t35 * t27;
t57 = t35 * t31;
t56 = rSges(4,3) - t31;
t19 = t24 * qJ(4);
t55 = t19 + t20;
t54 = qJ(4) * t25;
t53 = rSges(3,3) + qJ(2);
t52 = -m(5) + t72;
t30 = cos(pkin(9));
t23 = pkin(2) * t30 + pkin(1);
t16 = t35 * t23;
t51 = pkin(3) * t62 + t35 * t19 + t16;
t50 = -t23 - t20;
t49 = pkin(4) * t64 + qJ(5) * t65 + t55;
t32 = sin(qJ(6));
t34 = cos(qJ(6));
t7 = t25 * t60 + t61;
t8 = t25 * t59 - t58;
t48 = t32 * t8 - t34 * t7;
t47 = -t32 * t7 - t34 * t8;
t46 = rSges(4,1) * t25 - rSges(4,2) * t24;
t44 = t27 * t34 - t29 * t32;
t43 = t27 * t32 + t29 * t34;
t42 = rSges(3,1) * t30 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t40 = -t8 * pkin(4) - t7 * qJ(5) - t57;
t39 = t50 - t19;
t10 = t25 * t61 + t60;
t9 = t25 * t58 - t59;
t38 = t10 * pkin(4) + t9 * qJ(5) + t51;
t15 = t35 * t54;
t12 = t33 * t54;
t3 = t10 * t34 + t32 * t9;
t2 = -t10 * t32 + t34 * t9;
t1 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t35) + g(2) * (rSges(2,1) * t35 - t33 * rSges(2,2))) - m(3) * ((g(1) * t53 + g(2) * t42) * t35 + (-g(1) * t42 + g(2) * t53) * t33) - m(4) * (g(2) * t16 + (g(1) * t56 + g(2) * t46) * t35 + (g(1) * (-t23 - t46) + g(2) * t56) * t33) - m(5) * (g(1) * (-t8 * rSges(5,1) + t7 * rSges(5,2) - t57) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + rSges(5,3) * t66 + t51) + (g(1) * (-rSges(5,3) * t24 + t39) - t70) * t33) - m(6) * (g(1) * (-t8 * rSges(6,1) - t7 * rSges(6,3) + t40) + g(2) * (t10 * rSges(6,1) + rSges(6,2) * t66 + t9 * rSges(6,3) + t38) + (g(1) * (-rSges(6,2) * t24 + t39) - t70) * t33) - m(7) * (g(1) * (t47 * rSges(7,1) + t48 * rSges(7,2) - t8 * pkin(5) + t40) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t10 * pkin(5) + t67 * t66 + t38) + (-t70 + (t50 + (-qJ(4) - t67) * t24) * g(1)) * t33) (-m(3) - m(4) + t52) * (g(1) * t33 - g(2) * t35) -m(4) * (g(3) * t46 + t73 * (-rSges(4,1) * t24 - rSges(4,2) * t25)) - m(5) * (g(1) * (rSges(5,3) * t62 + t15) + g(2) * (rSges(5,3) * t63 + t12) + g(3) * (rSges(5,1) * t64 - rSges(5,2) * t65 + t55) + (g(3) * rSges(5,3) + t73 * (-rSges(5,1) * t29 + rSges(5,2) * t27 - pkin(3))) * t24) - m(6) * (g(1) * (rSges(6,2) * t62 + t15) + g(2) * (rSges(6,2) * t63 + t12) + g(3) * (rSges(6,1) * t64 + rSges(6,3) * t65 + t49) + (g(3) * rSges(6,2) + t73 * (-pkin(3) + (-rSges(6,1) - pkin(4)) * t29 + (-rSges(6,3) - qJ(5)) * t27)) * t24) - m(7) * (g(1) * t15 + g(2) * t12 + g(3) * t49 + (g(3) * (t43 * rSges(7,1) + t44 * rSges(7,2) + t29 * pkin(5)) + t73 * t67) * t25 + (g(3) * t67 + t73 * (-pkin(3) + (-t32 * rSges(7,1) - t34 * rSges(7,2) - qJ(5)) * t27 + (-t34 * rSges(7,1) + t32 * rSges(7,2) - pkin(4) - pkin(5)) * t29)) * t24) t52 * (-g(3) * t25 + t73 * t24) t72 * (g(1) * t9 + g(2) * t7 + t27 * t68) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (-t48 * rSges(7,1) + t47 * rSges(7,2)) + (t44 * rSges(7,1) - t43 * rSges(7,2)) * t68)];
taug  = t1(:);
