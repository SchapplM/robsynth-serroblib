% Calculate Gravitation load on the joints for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:08
% EndTime: 2019-03-09 02:33:09
% DurationCPUTime: 0.53s
% Computational Cost: add. (309->106), mult. (301->137), div. (0->0), fcn. (258->10), ass. (0->60)
t31 = pkin(10) + qJ(4);
t26 = qJ(5) + t31;
t23 = cos(t26);
t78 = rSges(7,3) + pkin(9);
t79 = t78 * t23;
t22 = sin(t26);
t77 = t78 * t22;
t38 = cos(qJ(1));
t69 = g(2) * t38;
t36 = sin(qJ(1));
t76 = g(1) * t38 + g(2) * t36;
t35 = sin(qJ(6));
t64 = rSges(7,2) * t35;
t50 = t23 * t64;
t75 = t50 * t69;
t24 = sin(t31);
t74 = pkin(4) * t24;
t25 = cos(t31);
t73 = pkin(4) * t25;
t72 = g(1) * t36;
t32 = sin(pkin(10));
t68 = t32 * pkin(3);
t66 = rSges(6,1) * t23;
t37 = cos(qJ(6));
t65 = rSges(7,1) * t37;
t63 = t22 * t36;
t62 = t22 * t38;
t61 = t23 * t36;
t60 = t36 * t35;
t59 = t36 * t37;
t58 = t38 * t35;
t57 = t38 * t37;
t34 = -pkin(7) - qJ(3);
t56 = rSges(5,3) - t34;
t55 = t38 * pkin(1) + t36 * qJ(2);
t54 = rSges(4,3) + qJ(3);
t28 = t38 * qJ(2);
t30 = -pkin(8) + t34;
t9 = t68 + t74;
t53 = t36 * t30 + t38 * t9 + t28;
t52 = t36 * t9 + t55;
t51 = t23 * t65;
t49 = -m(4) - m(5) - m(6) - m(7);
t48 = -pkin(5) - t65;
t47 = rSges(6,1) * t61 - rSges(6,2) * t63;
t46 = rSges(4,1) * t32 + rSges(4,2) * cos(pkin(10));
t45 = rSges(5,1) * t25 - rSges(5,2) * t24;
t44 = -t24 * rSges(5,1) - t25 * rSges(5,2);
t43 = t22 * rSges(6,1) + t23 * rSges(6,2);
t42 = g(1) * t28 + g(2) * t55;
t41 = -t44 + t68;
t40 = pkin(5) * t61 + t78 * t63 + (-t50 + t51) * t36;
t39 = t79 + (t48 + t64) * t22;
t18 = t36 * t73;
t13 = rSges(6,2) * t62;
t4 = t22 * t57 - t60;
t3 = t22 * t58 + t59;
t2 = t22 * t59 + t58;
t1 = -t22 * t60 + t57;
t5 = [-m(2) * (g(1) * (-t36 * rSges(2,1) - t38 * rSges(2,2)) + g(2) * (t38 * rSges(2,1) - t36 * rSges(2,2))) - m(3) * (g(1) * (t38 * rSges(3,3) + t28 + (rSges(3,2) - pkin(1)) * t36) + g(2) * (-t38 * rSges(3,2) + t36 * rSges(3,3) + t55)) - m(4) * ((g(1) * t46 + g(2) * t54) * t38 + (g(1) * (-pkin(1) - t54) + g(2) * t46) * t36 + t42) - m(5) * ((g(1) * t41 + g(2) * t56) * t38 + (g(1) * (-pkin(1) - t56) + g(2) * t41) * t36 + t42) - m(6) * (g(1) * t53 + g(2) * t52 + (g(1) * t43 + g(2) * (rSges(6,3) - t30)) * t38 + (g(1) * (-rSges(6,3) - pkin(1)) + g(2) * t43) * t36) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) - t36 * pkin(1) + pkin(5) * t62 + t53) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t63 - t38 * t30 + t52) - t76 * t79) (-m(3) + t49) * (-t69 + t72) t49 * t76, -m(5) * (g(3) * t44 + t45 * t72) - m(6) * (g(1) * (t18 + t47) + g(2) * t13 + g(3) * (-t43 - t74)) - m(7) * (g(1) * (t18 + t40) + t75 + g(3) * (t39 - t74)) + (m(5) * t45 - m(6) * (-t66 - t73) - m(7) * (-pkin(5) * t23 - t51 - t73 - t77)) * t69, -m(6) * (g(1) * t47 + g(2) * (-t38 * t66 + t13) - g(3) * t43) - m(7) * (g(1) * t40 + t75 + g(3) * t39 + (t48 * t23 - t77) * t69) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t4 * rSges(7,2)) + g(3) * (-rSges(7,1) * t35 - rSges(7,2) * t37) * t23)];
taug  = t5(:);
