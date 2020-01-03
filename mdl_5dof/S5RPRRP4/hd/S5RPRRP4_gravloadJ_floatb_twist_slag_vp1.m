% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:15
% EndTime: 2020-01-03 11:49:17
% DurationCPUTime: 0.45s
% Computational Cost: add. (233->106), mult. (324->142), div. (0->0), fcn. (316->8), ass. (0->47)
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t66 = g(2) * t38 + g(3) * t36;
t39 = -pkin(7) - pkin(6);
t34 = cos(pkin(8));
t32 = qJ(3) + qJ(4);
t25 = sin(t32);
t48 = t38 * t25;
t26 = cos(t32);
t51 = t36 * t26;
t10 = t34 * t51 - t48;
t47 = t38 * t26;
t52 = t36 * t25;
t9 = -t34 * t52 - t47;
t65 = t9 * rSges(6,1) - t10 * rSges(6,2);
t64 = t9 * rSges(5,1) - t10 * rSges(5,2);
t11 = t34 * t48 - t51;
t12 = t34 * t47 + t52;
t63 = t11 * rSges(6,1) + t12 * rSges(6,2);
t62 = t11 * rSges(5,1) + t12 * rSges(5,2);
t61 = pkin(4) * t25;
t33 = sin(pkin(8));
t60 = g(1) * t33;
t35 = sin(qJ(3));
t57 = t35 * pkin(3);
t55 = rSges(3,2) * t33;
t54 = t34 * t36;
t53 = t34 * t38;
t50 = t36 * t35;
t37 = cos(qJ(3));
t49 = t36 * t37;
t46 = t38 * t35;
t45 = t38 * t37;
t29 = t37 * pkin(3);
t19 = pkin(4) * t26 + t29;
t44 = t38 * pkin(1) + t36 * qJ(2);
t43 = -rSges(5,1) * t25 - rSges(5,2) * t26;
t42 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t15 = t34 * t46 - t49;
t13 = -t34 * t50 - t45;
t41 = (t29 + pkin(2)) * t34 + (rSges(5,3) - t39) * t33;
t40 = (pkin(2) + t19) * t34 + (rSges(6,3) + qJ(5) - t39) * t33;
t28 = t36 * pkin(1);
t18 = t57 + t61;
t16 = t34 * t45 + t50;
t14 = t34 * t49 - t46;
t1 = [-m(2) * (g(2) * (t38 * rSges(2,1) - t36 * rSges(2,2)) + g(3) * (t36 * rSges(2,1) + t38 * rSges(2,2))) - m(3) * (g(2) * (t36 * rSges(3,3) + t44) + g(3) * (rSges(3,1) * t54 - t36 * t55 + t28) + (g(2) * (rSges(3,1) * t34 - t55) + g(3) * (-rSges(3,3) - qJ(2))) * t38) - m(4) * (g(2) * (t16 * rSges(4,1) - t15 * rSges(4,2) + pkin(2) * t53 + t44) + g(3) * (t14 * rSges(4,1) + t13 * rSges(4,2) + pkin(2) * t54 - t38 * qJ(2) + t28) + t66 * t33 * (rSges(4,3) + pkin(6))) - m(5) * (g(2) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t44) + g(3) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t28) + (g(2) * t57 + g(3) * t41) * t36 + (g(2) * t41 + g(3) * (-qJ(2) - t57)) * t38) - m(6) * (g(2) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t44) + g(3) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t28) + (g(2) * t18 + g(3) * t40) * t36 + (g(2) * t40 + g(3) * (-qJ(2) - t18)) * t38), -(-m(3) - m(4) - m(5) - m(6)) * t66, -m(4) * (g(2) * (t13 * rSges(4,1) - t14 * rSges(4,2)) + g(3) * (t15 * rSges(4,1) + t16 * rSges(4,2))) - m(5) * (g(2) * (t13 * pkin(3) + t64) + g(3) * (t15 * pkin(3) + t62)) - m(6) * (g(2) * (-t18 * t54 - t38 * t19 + t65) + g(3) * (t18 * t53 - t36 * t19 + t63)) + (-m(4) * (-rSges(4,1) * t35 - rSges(4,2) * t37) - m(5) * (t43 - t57) - m(6) * (-t18 + t42)) * t60, -m(5) * (g(2) * t64 + g(3) * t62) - m(6) * (g(2) * (t9 * pkin(4) + t65) + g(3) * (t11 * pkin(4) + t63)) + (-m(5) * t43 - m(6) * (t42 - t61)) * t60, -m(6) * (-g(1) * t34 + (g(2) * t36 - g(3) * t38) * t33)];
taug = t1(:);
