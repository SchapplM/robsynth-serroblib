% Calculate Gravitation load on the joints for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:29
% EndTime: 2019-12-31 19:12:31
% DurationCPUTime: 0.48s
% Computational Cost: add. (194->95), mult. (259->126), div. (0->0), fcn. (224->8), ass. (0->49)
t25 = qJ(3) + qJ(4);
t21 = cos(t25);
t63 = rSges(6,3) + pkin(8);
t64 = t63 * t21;
t20 = sin(t25);
t62 = t63 * t20;
t31 = cos(qJ(1));
t58 = g(2) * t31;
t26 = sin(qJ(5));
t52 = rSges(6,2) * t26;
t42 = t21 * t52;
t61 = t42 * t58;
t30 = cos(qJ(3));
t60 = pkin(3) * t30;
t28 = sin(qJ(1));
t59 = g(1) * t28;
t27 = sin(qJ(3));
t57 = t27 * pkin(3);
t56 = rSges(4,3) + pkin(6);
t54 = rSges(5,1) * t21;
t29 = cos(qJ(5));
t53 = rSges(6,1) * t29;
t51 = t20 * t28;
t50 = t20 * t31;
t49 = t21 * t28;
t48 = t28 * t26;
t47 = t28 * t29;
t46 = t31 * t26;
t45 = t31 * t29;
t44 = t31 * pkin(1) + t28 * qJ(2);
t43 = t21 * t53;
t41 = t28 * t57 + t44;
t23 = t31 * qJ(2);
t32 = -pkin(7) - pkin(6);
t40 = t28 * t32 + t31 * t57 + t23;
t39 = -pkin(4) - t53;
t38 = rSges(5,1) * t49 - rSges(5,2) * t51;
t37 = rSges(4,1) * t30 - rSges(4,2) * t27;
t36 = t27 * rSges(4,1) + t30 * rSges(4,2);
t35 = t20 * rSges(5,1) + t21 * rSges(5,2);
t34 = pkin(4) * t49 + t63 * t51 + (-t42 + t43) * t28;
t33 = t64 + (t39 + t52) * t20;
t15 = t28 * t60;
t11 = rSges(5,2) * t50;
t4 = t20 * t45 - t48;
t3 = t20 * t46 + t47;
t2 = t20 * t47 + t46;
t1 = -t20 * t48 + t45;
t5 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - t31 * rSges(2,2)) + g(2) * (t31 * rSges(2,1) - t28 * rSges(2,2))) - m(3) * (g(1) * (t31 * rSges(3,3) + t23 + (rSges(3,2) - pkin(1)) * t28) + g(2) * (-t31 * rSges(3,2) + t28 * rSges(3,3) + t44)) - m(4) * (g(1) * t23 + g(2) * t44 + (g(1) * t36 + g(2) * t56) * t31 + (g(1) * (-pkin(1) - t56) + g(2) * t36) * t28) - m(5) * (g(1) * t40 + g(2) * t41 + (g(1) * t35 + g(2) * (rSges(5,3) - t32)) * t31 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t35) * t28) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t28 * pkin(1) + pkin(4) * t50 + t40) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t51 - t31 * t32 + t41) - (g(1) * t31 + g(2) * t28) * t64), (-m(3) - m(4) - m(5) - m(6)) * (-t58 + t59), -m(4) * (-g(3) * t36 + t37 * t59) - m(5) * (g(1) * (t15 + t38) + g(2) * t11 + g(3) * (-t35 - t57)) - m(6) * (g(1) * (t15 + t34) + t61 + g(3) * (t33 - t57)) + (m(4) * t37 - m(5) * (-t54 - t60) - m(6) * (-pkin(4) * t21 - t43 - t60 - t62)) * t58, -m(5) * (g(1) * t38 + g(2) * (-t31 * t54 + t11) - g(3) * t35) - m(6) * (g(1) * t34 + t61 + g(3) * t33 + (t39 * t21 - t62) * t58), -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t4 * rSges(6,2)) + g(3) * (-rSges(6,1) * t26 - rSges(6,2) * t29) * t21)];
taug = t5(:);
