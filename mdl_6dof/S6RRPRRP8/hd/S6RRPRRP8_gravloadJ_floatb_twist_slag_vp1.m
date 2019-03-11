% Calculate Gravitation load on the joints for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:00
% EndTime: 2019-03-09 12:22:02
% DurationCPUTime: 0.89s
% Computational Cost: add. (580->152), mult. (610->208), div. (0->0), fcn. (597->10), ass. (0->66)
t84 = rSges(7,1) + pkin(5);
t43 = pkin(10) + qJ(4);
t35 = sin(t43);
t36 = cos(t43);
t48 = sin(qJ(1));
t49 = cos(qJ(2));
t50 = cos(qJ(1));
t78 = t50 * t49;
t17 = -t35 * t78 + t48 * t36;
t71 = rSges(7,3) + qJ(6);
t46 = -pkin(8) - qJ(3);
t76 = rSges(5,3) - t46;
t87 = g(2) * t48;
t95 = g(1) * t50 + t87;
t37 = qJ(5) + t43;
t32 = sin(t37);
t33 = cos(t37);
t94 = t71 * t32 + t84 * t33;
t79 = t48 * t49;
t11 = t32 * t79 + t33 * t50;
t12 = -t50 * t32 + t33 * t79;
t93 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t32 * t78 - t48 * t33;
t14 = t48 * t32 + t33 * t78;
t92 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t44 = sin(pkin(10));
t91 = pkin(3) * t44;
t90 = pkin(4) * t35;
t89 = g(1) * t48;
t45 = cos(pkin(10));
t34 = t45 * pkin(3) + pkin(2);
t24 = pkin(4) * t36 + t34;
t21 = t49 * t24;
t86 = g(3) * t21;
t47 = sin(qJ(2));
t85 = g(3) * t47;
t83 = rSges(3,2) * t47;
t82 = t32 * t47;
t42 = -pkin(9) + t46;
t77 = rSges(7,2) - t42;
t75 = rSges(6,3) - t42;
t74 = t71 * t33 * t47;
t73 = t50 * pkin(1) + t48 * pkin(7);
t72 = rSges(4,3) + qJ(3);
t25 = t90 + t91;
t40 = t50 * pkin(7);
t69 = t48 * t47 * t42 + t50 * t25 + t40;
t68 = -pkin(1) - t21;
t67 = t77 * t50;
t66 = t75 * t50;
t65 = t50 * t72;
t64 = t24 * t78 + t48 * t25 + t73;
t63 = rSges(3,1) * t49 - t83;
t61 = rSges(6,1) * t33 - rSges(6,2) * t32;
t60 = -rSges(6,1) * t32 - rSges(6,2) * t33;
t59 = -t84 * t11 + t71 * t12;
t58 = t17 * pkin(4);
t57 = -t84 * t13 + t71 * t14;
t56 = rSges(4,1) * t45 - rSges(4,2) * t44 + pkin(2);
t15 = t35 * t79 + t36 * t50;
t55 = rSges(5,1) * t36 - rSges(5,2) * t35 + t34;
t53 = t15 * pkin(4);
t52 = t34 * t49 + t47 * t76;
t18 = t48 * t35 + t36 * t78;
t16 = t35 * t50 - t36 * t79;
t1 = [-m(2) * (g(1) * (-t48 * rSges(2,1) - rSges(2,2) * t50) + g(2) * (rSges(2,1) * t50 - t48 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t50 + t40) + g(2) * (rSges(3,1) * t78 - t50 * t83 + t73) + (g(1) * (-pkin(1) - t63) + g(2) * rSges(3,3)) * t48) - m(4) * (g(1) * (-pkin(2) * t79 - t48 * pkin(1) + t40 + (t44 * t50 - t45 * t79) * rSges(4,1) + (t44 * t79 + t45 * t50) * rSges(4,2)) + g(2) * (pkin(2) * t78 + (t48 * t44 + t45 * t78) * rSges(4,1) + (-t44 * t78 + t48 * t45) * rSges(4,2) + t73) + (g(2) * t65 - t72 * t89) * t47) - m(5) * (g(1) * (t16 * rSges(5,1) + t15 * rSges(5,2) + t40) + g(2) * (t18 * rSges(5,1) + t17 * rSges(5,2) + t73) + (g(1) * t91 + g(2) * t52) * t50 + (g(1) * (-pkin(1) - t52) + g(2) * t91) * t48) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + t69) + g(2) * (t14 * rSges(6,1) - t13 * rSges(6,2) + t47 * t66 + t64) + (-rSges(6,3) * t47 + t68) * t89) - m(7) * (g(1) * (-t71 * t11 - t12 * t84 + t69) + g(2) * (t71 * t13 + t84 * t14 + t47 * t67 + t64) + (-rSges(7,2) * t47 + t68) * t89) -m(3) * (g(3) * t63 + t95 * (-rSges(3,1) * t47 - rSges(3,2) * t49)) - m(4) * ((g(1) * t65 + g(3) * t56 + t72 * t87) * t49 + (g(3) * t72 - t56 * t95) * t47) - m(5) * ((g(3) * t55 + t76 * t95) * t49 + (g(3) * t76 - t55 * t95) * t47) - m(6) * (t86 + (g(1) * t66 + g(3) * t61 + t75 * t87) * t49 + (g(3) * t75 + t95 * (-t24 - t61)) * t47) - m(7) * (t86 + (g(1) * t67 + g(3) * t94 + t77 * t87) * t49 + (g(3) * t77 + t95 * (-t24 - t94)) * t47) (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t49 + t47 * t95) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (t58 + t92) + g(2) * (-t53 + t93)) - m(7) * (g(1) * (t57 + t58) + g(2) * (-t53 + t59) + g(3) * t74) + (-m(5) * (-rSges(5,1) * t35 - rSges(5,2) * t36) - m(6) * (t60 - t90) - m(7) * (-t32 * t84 - t90)) * t85, -m(6) * (g(1) * t92 + g(2) * t93 + t60 * t85) - m(7) * (g(1) * t57 + g(2) * t59 + g(3) * (-t84 * t82 + t74)) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t82)];
taug  = t1(:);
