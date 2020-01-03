% Calculate Gravitation load on the joints for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:09
% DurationCPUTime: 0.56s
% Computational Cost: add. (279->102), mult. (334->138), div. (0->0), fcn. (293->8), ass. (0->60)
t32 = sin(qJ(5));
t35 = cos(qJ(5));
t87 = rSges(6,1) * t32 + rSges(6,2) * t35;
t31 = qJ(2) + qJ(3);
t29 = cos(t31);
t77 = rSges(6,3) + pkin(8);
t86 = t77 * t29;
t85 = t87 * t29;
t28 = sin(t31);
t84 = t29 * rSges(4,1) - t28 * rSges(4,2);
t44 = -t29 * rSges(5,2) + t28 * rSges(5,3);
t37 = cos(qJ(1));
t34 = sin(qJ(1));
t80 = g(2) * t34;
t83 = g(1) * t37 + t80;
t33 = sin(qJ(2));
t82 = pkin(2) * t33;
t79 = g(3) * t29;
t26 = t29 * pkin(3);
t78 = rSges(3,3) + pkin(6);
t38 = -pkin(7) - pkin(6);
t76 = pkin(4) - t38;
t72 = t28 * t34;
t71 = t28 * t37;
t69 = t29 * t37;
t68 = t34 * t32;
t67 = t34 * t35;
t66 = t37 * t32;
t65 = t37 * t35;
t64 = rSges(5,1) - t38;
t63 = rSges(4,3) - t38;
t21 = t28 * qJ(4);
t62 = t21 + t26;
t61 = qJ(4) * t29;
t60 = -pkin(3) - t77;
t10 = t34 * t61;
t59 = t85 * t34 + t10;
t12 = t37 * t61;
t58 = t85 * t37 + t12;
t55 = t34 * t29 * rSges(5,3) + rSges(5,2) * t72 + t10;
t36 = cos(qJ(2));
t30 = t36 * pkin(2);
t27 = t30 + pkin(1);
t13 = t37 * t27;
t54 = pkin(3) * t69 + t37 * t21 + t13;
t53 = rSges(5,2) * t71 + rSges(5,3) * t69 + t12;
t51 = -t27 - t21;
t50 = g(1) * t60;
t49 = -pkin(3) * t28 - t82;
t48 = t36 * rSges(3,1) - t33 * rSges(3,2);
t45 = -rSges(4,1) * t28 - rSges(4,2) * t29;
t43 = t62 + t44;
t42 = t87 * t28 + t62 + t86;
t41 = pkin(1) + t48;
t39 = (t37 * t50 + t60 * t80) * t28;
t5 = -t28 * t68 + t65;
t4 = t28 * t67 + t66;
t3 = t28 * t66 + t67;
t2 = t28 * t65 - t68;
t1 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) - t34 * rSges(2,2))) - m(3) * ((g(1) * t78 + g(2) * t41) * t37 + (-g(1) * t41 + g(2) * t78) * t34) - m(4) * (g(2) * t13 + (g(1) * t63 + g(2) * t84) * t37 + (g(1) * (-t27 - t84) + g(2) * t63) * t34) - m(5) * (g(2) * t54 + (g(1) * t64 + g(2) * t44) * t37 + (g(1) * (-t44 + t51 - t26) + g(2) * t64) * t34) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t54) + (g(1) * t76 + g(2) * t86) * t37 + (g(1) * t51 + g(2) * t76 + t29 * t50) * t34), -m(3) * (g(3) * t48 + t83 * (-rSges(3,1) * t33 - rSges(3,2) * t36)) - m(4) * (g(3) * (t30 + t84) + t83 * (t45 - t82)) - m(5) * (g(1) * (t37 * t49 + t53) + g(2) * (t34 * t49 + t55) + g(3) * (t30 + t43)) - m(6) * (g(1) * (-t37 * t82 + t58) + g(2) * (-t34 * t82 + t59) + g(3) * (t30 + t42) + t39), -m(4) * (g(3) * t84 + t83 * t45) - m(5) * (g(1) * (-pkin(3) * t71 + t53) + g(2) * (-pkin(3) * t72 + t55) + g(3) * t43) - m(6) * (g(1) * t58 + g(2) * t59 + g(3) * t42 + t39), (-m(5) - m(6)) * (t83 * t28 - t79), -m(6) * (g(1) * (t2 * rSges(6,1) - t3 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t5 * rSges(6,2)) + (-rSges(6,1) * t35 + rSges(6,2) * t32) * t79)];
taug = t1(:);
