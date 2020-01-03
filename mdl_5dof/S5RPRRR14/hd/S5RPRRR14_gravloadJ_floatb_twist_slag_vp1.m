% Calculate Gravitation load on the joints for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:54
% EndTime: 2019-12-31 19:16:58
% DurationCPUTime: 1.24s
% Computational Cost: add. (699->140), mult. (1875->217), div. (0->0), fcn. (2380->14), ass. (0->64)
t52 = cos(qJ(1));
t86 = cos(pkin(11));
t88 = cos(pkin(5));
t77 = t88 * t86;
t83 = sin(pkin(11));
t92 = sin(qJ(1));
t66 = -t52 * t77 + t92 * t83;
t84 = sin(pkin(6));
t85 = sin(pkin(5));
t74 = t85 * t84;
t87 = cos(pkin(6));
t101 = t52 * t74 + t66 * t87;
t76 = t88 * t83;
t35 = t52 * t76 + t92 * t86;
t49 = sin(qJ(3));
t93 = cos(qJ(3));
t15 = t101 * t93 + t35 * t49;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t18 = t101 * t49 - t35 * t93;
t75 = t85 * t87;
t27 = -t52 * t75 + t66 * t84;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t6 = t18 * t51 - t27 * t48;
t105 = t15 * t50 + t6 * t47;
t104 = -t15 * t47 + t6 * t50;
t100 = t18 * t48 + t27 * t51;
t60 = t52 * t83 + t92 * t77;
t53 = t60 * t84 + t92 * t75;
t98 = t60 * t87 - t92 * t74;
t97 = t86 * t75 + t84 * t88;
t96 = t51 * pkin(4);
t95 = rSges(6,3) + pkin(10);
t94 = pkin(9) + rSges(5,3);
t91 = t47 * t51;
t90 = t50 * t51;
t79 = t85 * t92;
t89 = t52 * pkin(1) + qJ(2) * t79;
t82 = t52 * t85;
t81 = -t92 * pkin(1) + qJ(2) * t82;
t78 = -rSges(5,1) * t51 + rSges(5,2) * t48;
t73 = t85 * t83;
t71 = rSges(6,1) * t50 - rSges(6,2) * t47 + pkin(4);
t57 = -t35 * pkin(2) - t27 * pkin(8) + t81;
t56 = t18 * pkin(3) + t57;
t36 = t52 * t86 - t92 * t76;
t55 = t36 * pkin(2) + t53 * pkin(8) + t89;
t20 = t36 * t93 - t49 * t98;
t54 = t20 * pkin(3) + t55;
t34 = -t86 * t74 + t88 * t87;
t25 = t49 * t97 + t93 * t73;
t24 = t49 * t73 - t93 * t97;
t23 = t24 * pkin(3);
t19 = t36 * t49 + t93 * t98;
t14 = t25 * t51 + t34 * t48;
t13 = -t25 * t48 + t34 * t51;
t11 = t19 * pkin(3);
t9 = t15 * pkin(3);
t8 = t20 * t51 + t53 * t48;
t7 = t20 * t48 - t53 * t51;
t2 = t19 * t47 + t8 * t50;
t1 = t19 * t50 - t8 * t47;
t3 = [-m(2) * (g(1) * (-t92 * rSges(2,1) - t52 * rSges(2,2)) + g(2) * (t52 * rSges(2,1) - t92 * rSges(2,2))) - m(3) * (g(1) * (-t35 * rSges(3,1) + t66 * rSges(3,2) + rSges(3,3) * t82 + t81) + g(2) * (t36 * rSges(3,1) - t60 * rSges(3,2) + rSges(3,3) * t79 + t89)) - m(4) * (g(1) * (t18 * rSges(4,1) + rSges(4,2) * t15 - rSges(4,3) * t27 + t57) + g(2) * (t20 * rSges(4,1) - t19 * rSges(4,2) + t53 * rSges(4,3) + t55)) - m(5) * (g(1) * (t6 * rSges(5,1) - rSges(5,2) * t100 - t15 * t94 + t56) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t94 * t19 + t54)) - m(6) * (g(1) * (t104 * rSges(6,1) - rSges(6,2) * t105 + t6 * pkin(4) - t15 * pkin(9) + t95 * t100 + t56) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) + t19 * pkin(9) + t95 * t7 + t54)), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t79 - g(2) * t82 + g(3) * t88), -m(4) * (g(1) * (-t19 * rSges(4,1) - t20 * rSges(4,2)) + g(2) * (-t15 * rSges(4,1) + rSges(4,2) * t18) + g(3) * (-t24 * rSges(4,1) - t25 * rSges(4,2))) - m(5) * (g(1) * (t78 * t19 + t94 * t20 - t11) + g(2) * (t78 * t15 - t18 * t94 - t9) + g(3) * (t78 * t24 + t94 * t25 - t23)) + (-g(1) * (-t19 * t96 - t11 + t20 * pkin(9) + (-t19 * t90 + t20 * t47) * rSges(6,1) + (t19 * t91 + t20 * t50) * rSges(6,2)) - g(2) * (-t15 * t96 - t9 - t18 * pkin(9) + (-t15 * t90 - t18 * t47) * rSges(6,1) + (t15 * t91 - t18 * t50) * rSges(6,2)) - g(3) * (-t24 * t96 - t23 + t25 * pkin(9) + (-t24 * t90 + t25 * t47) * rSges(6,1) + (t24 * t91 + t25 * t50) * rSges(6,2)) - (-g(1) * t19 - g(2) * t15 - g(3) * t24) * t48 * t95) * m(6), -m(5) * (g(1) * (-t7 * rSges(5,1) - t8 * rSges(5,2)) + g(2) * (rSges(5,1) * t100 + rSges(5,2) * t6) + g(3) * (t13 * rSges(5,1) - t14 * rSges(5,2))) - m(6) * (g(1) * (-t71 * t7 + t95 * t8) + (t71 * t13 + t95 * t14) * g(3) + (t100 * t71 - t6 * t95) * g(2)), -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (rSges(6,1) * t105 + t104 * rSges(6,2)) + g(3) * ((-t14 * t47 + t24 * t50) * rSges(6,1) + (-t14 * t50 - t24 * t47) * rSges(6,2)))];
taug = t3(:);
