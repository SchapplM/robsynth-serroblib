% Calculate Gravitation load on the joints for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:50:56
% DurationCPUTime: 0.27s
% Computational Cost: add. (268->70), mult. (208->85), div. (0->0), fcn. (169->8), ass. (0->35)
t57 = rSges(6,1) + pkin(4);
t27 = pkin(8) + qJ(4);
t23 = cos(t27);
t61 = t57 * t23;
t22 = sin(t27);
t60 = rSges(5,1) * t23 - t22 * rSges(5,2);
t59 = qJ(3) + rSges(4,3);
t30 = cos(pkin(8));
t58 = rSges(4,2) * sin(pkin(8)) - pkin(2) - rSges(4,1) * t30;
t28 = qJ(1) + qJ(2);
t24 = sin(t28);
t25 = cos(t28);
t56 = g(1) * t25 + g(2) * t24;
t55 = rSges(6,3) + qJ(5);
t32 = sin(qJ(1));
t54 = pkin(1) * t32;
t47 = t22 * t25;
t46 = t23 * t25;
t31 = -pkin(7) - qJ(3);
t45 = t25 * t31;
t44 = qJ(5) * t22;
t43 = t25 * rSges(3,1) - rSges(3,2) * t24;
t42 = t25 * rSges(6,2) - t45;
t20 = pkin(3) * t30 + pkin(2);
t3 = t25 * t20;
t41 = t24 * rSges(6,2) + rSges(6,3) * t47 + t25 * t44 + t57 * t46 + t3;
t40 = -rSges(3,1) * t24 - rSges(3,2) * t25;
t39 = t59 * t24 - t58 * t25;
t38 = t58 * t24 + t59 * t25;
t37 = rSges(5,1) * t46 - rSges(5,2) * t47 + t3 + (rSges(5,3) - t31) * t24;
t36 = t25 * rSges(5,3) - t45 + (-t20 - t60) * t24;
t35 = (g(1) * (-rSges(6,3) * t22 - t20 - t44 - t61) - g(2) * t31) * t24;
t33 = cos(qJ(1));
t26 = t33 * pkin(1);
t1 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - rSges(2,2) * t33) + g(2) * (rSges(2,1) * t33 - t32 * rSges(2,2))) - m(3) * (g(1) * (t40 - t54) + g(2) * (t26 + t43)) - m(4) * (g(1) * (t38 - t54) + g(2) * (t26 + t39)) - m(5) * (g(1) * (t36 - t54) + g(2) * (t26 + t37)) - m(6) * (g(1) * (t42 - t54) + g(2) * (t26 + t41) + t35), -m(3) * (g(1) * t40 + g(2) * t43) - m(4) * (g(1) * t38 + g(2) * t39) - m(5) * (g(1) * t36 + g(2) * t37) - m(6) * (g(1) * t42 + g(2) * t41 + t35), (-m(4) - m(5) - m(6)) * (g(1) * t24 - g(2) * t25), (-m(5) * t60 - m(6) * (t55 * t22 + t61)) * g(3) + t56 * (-m(5) * (-rSges(5,1) * t22 - rSges(5,2) * t23) - m(6) * (-t57 * t22 + t55 * t23)), -m(6) * (-g(3) * t23 + t56 * t22)];
taug = t1(:);
