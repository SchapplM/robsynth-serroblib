% Calculate Gravitation load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:53
% EndTime: 2019-12-31 20:16:55
% DurationCPUTime: 0.45s
% Computational Cost: add. (297->89), mult. (283->120), div. (0->0), fcn. (245->10), ass. (0->52)
t78 = rSges(6,3) + pkin(8);
t31 = qJ(2) + pkin(9);
t28 = qJ(4) + t31;
t23 = sin(t28);
t24 = cos(t28);
t33 = sin(qJ(5));
t67 = rSges(6,2) * t33;
t77 = t23 * t67 + t24 * t78;
t75 = t24 * rSges(5,1) - t23 * rSges(5,2);
t26 = sin(t31);
t27 = cos(t31);
t37 = cos(qJ(2));
t29 = t37 * pkin(2);
t74 = t27 * rSges(4,1) - t26 * rSges(4,2) + t29;
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t73 = g(1) * t38 + g(2) * t35;
t42 = t24 * pkin(4) + t78 * t23;
t34 = sin(qJ(2));
t70 = t34 * pkin(2);
t69 = rSges(3,3) + pkin(6);
t36 = cos(qJ(5));
t68 = rSges(6,1) * t36;
t63 = t35 * t33;
t62 = t35 * t36;
t61 = t38 * t33;
t60 = t38 * t36;
t32 = -qJ(3) - pkin(6);
t59 = rSges(4,3) - t32;
t30 = -pkin(7) + t32;
t58 = rSges(5,3) - t30;
t57 = pkin(3) * t27 + t29;
t56 = t77 * t35;
t55 = t77 * t38;
t51 = t37 * rSges(3,1) - t34 * rSges(3,2);
t47 = -rSges(5,1) * t23 - rSges(5,2) * t24;
t46 = pkin(1) + t51;
t45 = pkin(1) + t74;
t44 = t47 * t35;
t43 = t47 * t38;
t41 = t42 + (-t67 + t68) * t24;
t39 = t73 * (-pkin(4) - t68) * t23;
t11 = -pkin(3) * t26 - t70;
t10 = pkin(1) + t57;
t7 = t38 * t11;
t6 = t35 * t11;
t5 = t38 * t10;
t4 = t24 * t60 + t63;
t3 = -t24 * t61 + t62;
t2 = -t24 * t62 + t61;
t1 = t24 * t63 + t60;
t8 = [-m(2) * (g(1) * (-t35 * rSges(2,1) - t38 * rSges(2,2)) + g(2) * (t38 * rSges(2,1) - t35 * rSges(2,2))) - m(3) * ((g(1) * t69 + g(2) * t46) * t38 + (-g(1) * t46 + g(2) * t69) * t35) - m(4) * ((g(1) * t59 + g(2) * t45) * t38 + (-g(1) * t45 + g(2) * t59) * t35) - m(5) * (g(2) * t5 + (g(1) * t58 + g(2) * t75) * t38 + (g(1) * (-t10 - t75) + g(2) * t58) * t35) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t30 + g(2) * t42) * t38 + (g(1) * (-t10 - t42) - g(2) * t30) * t35), -m(3) * (g(3) * t51 + t73 * (-rSges(3,1) * t34 - rSges(3,2) * t37)) - m(4) * (g(3) * t74 + t73 * (-rSges(4,1) * t26 - rSges(4,2) * t27 - t70)) - m(5) * (g(1) * (t7 + t43) + g(2) * (t6 + t44) + g(3) * (t75 + t57)) - m(6) * (g(1) * (t7 + t55) + g(2) * (t6 + t56) + g(3) * (t41 + t57) + t39), (-m(4) - m(5) - m(6)) * (g(1) * t35 - g(2) * t38), -m(5) * (g(1) * t43 + g(2) * t44 + g(3) * t75) - m(6) * (g(1) * t55 + g(2) * t56 + g(3) * t41 + t39), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t33 - rSges(6,2) * t36) * t23)];
taug = t8(:);
