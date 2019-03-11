% Calculate Gravitation load on the joints for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:25
% EndTime: 2019-03-08 20:22:27
% DurationCPUTime: 0.88s
% Computational Cost: add. (662->134), mult. (1354->211), div. (0->0), fcn. (1671->14), ass. (0->69)
t127 = rSges(7,3) + pkin(10);
t54 = sin(pkin(12));
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t89 = cos(pkin(12));
t40 = -t64 * t54 - t61 * t89;
t55 = sin(pkin(11));
t57 = cos(pkin(11));
t58 = cos(pkin(6));
t96 = t58 * t64;
t75 = -t55 * t96 - t57 * t61;
t121 = t75 * pkin(2);
t73 = -t61 * t54 + t64 * t89;
t68 = t73 * t58;
t23 = t57 * t40 - t55 * t68;
t90 = t40 * t58;
t24 = t55 * t90 + t57 * t73;
t63 = cos(qJ(4));
t50 = t63 * pkin(4) + pkin(3);
t65 = -pkin(9) - pkin(8);
t126 = t23 * t50 - t24 * t65 + t121;
t125 = -t55 * t61 + t57 * t96;
t56 = sin(pkin(6));
t101 = t56 * t63;
t60 = sin(qJ(4));
t124 = t55 * t101 - t24 * t60;
t38 = t40 * t56;
t123 = t38 * t60 + t58 * t63;
t59 = sin(qJ(6));
t62 = cos(qJ(6));
t122 = t62 * rSges(7,1) - rSges(7,2) * t59 + pkin(5);
t120 = rSges(5,3) + pkin(8);
t20 = t55 * t40 + t57 * t68;
t37 = t73 * t56;
t119 = g(1) * t23 + g(2) * t20 + g(3) * t37;
t19 = -t55 * t73 + t57 * t90;
t103 = t56 * t57;
t53 = qJ(4) + qJ(5);
t51 = sin(t53);
t52 = cos(t53);
t11 = -t103 * t52 + t19 * t51;
t12 = -t103 * t51 - t19 * t52;
t118 = t11 * rSges(6,1) - t12 * rSges(6,2);
t105 = t55 * t56;
t13 = t105 * t52 - t24 * t51;
t14 = t105 * t51 + t24 * t52;
t117 = t13 * rSges(6,1) - t14 * rSges(6,2);
t113 = t52 * pkin(5);
t108 = t52 * t59;
t107 = t52 * t62;
t102 = t56 * t60;
t98 = t58 * t61;
t30 = t38 * t51 + t58 * t52;
t31 = -t38 * t52 + t58 * t51;
t91 = t30 * rSges(6,1) - t31 * rSges(6,2);
t47 = t56 * t64 * pkin(2);
t85 = t37 * t50 + t38 * t65 + t47;
t84 = -m(4) - m(5) - m(6) - m(7);
t81 = t125 * pkin(2);
t80 = t124 * pkin(4);
t79 = t123 * pkin(4);
t78 = rSges(6,1) * t52 - rSges(6,2) * t51;
t76 = -t101 * t57 + t19 * t60;
t74 = t19 * t65 + t20 * t50 + t81;
t72 = t76 * pkin(4);
t69 = t122 * t11 + t127 * t12;
t67 = t122 * t13 + t127 * t14;
t66 = t122 * t30 + t127 * t31;
t1 = [(-m(2) - m(3) + t84) * g(3), -m(3) * (g(1) * (t75 * rSges(3,1) + (t55 * t98 - t57 * t64) * rSges(3,2)) + g(2) * (t125 * rSges(3,1) + (-t55 * t64 - t57 * t98) * rSges(3,2)) + g(3) * (rSges(3,1) * t64 - rSges(3,2) * t61) * t56) - m(4) * (g(1) * (t23 * rSges(4,1) - rSges(4,2) * t24 + t121) + g(2) * (t20 * rSges(4,1) + t19 * rSges(4,2) + t81) + g(3) * (t37 * rSges(4,1) + t38 * rSges(4,2) + t47)) - m(5) * (g(1) * (t120 * t24 + t121) + g(2) * (-t120 * t19 + t81) + g(3) * (-t120 * t38 + t47) + t119 * (t63 * rSges(5,1) - t60 * rSges(5,2) + pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t24 + t23 * t78 + t126) + g(2) * (-t19 * rSges(6,3) + t20 * t78 + t74) + g(3) * (-t38 * rSges(6,3) + t37 * t78 + t85)) - m(7) * (g(1) * (t23 * t113 + (t107 * t23 + t24 * t59) * rSges(7,1) + (-t108 * t23 + t24 * t62) * rSges(7,2) + t126) + g(2) * (t20 * t113 + (t107 * t20 - t19 * t59) * rSges(7,1) + (-t108 * t20 - t19 * t62) * rSges(7,2) + t74) + g(3) * (t37 * t113 + (t107 * t37 - t38 * t59) * rSges(7,1) + (-t108 * t37 - t38 * t62) * rSges(7,2) + t85) + t119 * t51 * t127) t84 * (g(3) * t58 + (g(1) * t55 - g(2) * t57) * t56) -m(5) * (g(1) * (t124 * rSges(5,1) + (-t102 * t55 - t24 * t63) * rSges(5,2)) + g(2) * (t76 * rSges(5,1) + (t102 * t57 + t19 * t63) * rSges(5,2)) + g(3) * (t123 * rSges(5,1) + (t38 * t63 - t58 * t60) * rSges(5,2))) - m(6) * (g(1) * (t80 + t117) + g(2) * (t72 + t118) + g(3) * (t79 + t91)) - m(7) * (g(1) * (t67 + t80) + g(2) * (t72 + t69) + g(3) * (t66 + t79)) -m(6) * (g(1) * t117 + g(2) * t118 + g(3) * t91) - m(7) * (g(1) * t67 + g(2) * t69 + g(3) * t66) -m(7) * (g(1) * ((-t14 * t59 - t23 * t62) * rSges(7,1) + (-t14 * t62 + t23 * t59) * rSges(7,2)) + g(2) * ((-t12 * t59 - t20 * t62) * rSges(7,1) + (-t12 * t62 + t20 * t59) * rSges(7,2)) + g(3) * ((-t31 * t59 - t37 * t62) * rSges(7,1) + (-t31 * t62 + t37 * t59) * rSges(7,2)))];
taug  = t1(:);
