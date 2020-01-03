% Calculate Gravitation load on the joints for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:23
% EndTime: 2019-12-31 20:12:25
% DurationCPUTime: 0.50s
% Computational Cost: add. (180->101), mult. (378->140), div. (0->0), fcn. (361->6), ass. (0->42)
t49 = rSges(6,1) + pkin(4);
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t55 = g(1) * t25 + g(2) * t22;
t36 = rSges(6,3) + qJ(5);
t54 = -pkin(2) - pkin(7);
t53 = g(1) * t22;
t24 = cos(qJ(2));
t50 = g(3) * t24;
t16 = t24 * pkin(2);
t21 = sin(qJ(2));
t47 = t21 * t25;
t20 = sin(qJ(4));
t46 = t22 * t20;
t23 = cos(qJ(4));
t45 = t22 * t23;
t44 = t24 * rSges(4,2);
t43 = t24 * t25;
t42 = t25 * t20;
t41 = t25 * t23;
t12 = t21 * qJ(3);
t40 = t12 + t16;
t39 = t25 * pkin(1) + t22 * pkin(6);
t17 = t25 * pkin(6);
t38 = t25 * pkin(3) + t17;
t37 = qJ(3) * t24;
t35 = -rSges(6,2) + t54;
t34 = -rSges(5,3) + t54;
t33 = -pkin(1) - t12;
t32 = pkin(2) * t43 + t25 * t12 + t39;
t31 = t24 * rSges(3,1) - t21 * rSges(3,2);
t29 = rSges(5,1) * t20 + rSges(5,2) * t23;
t28 = t22 * pkin(3) + pkin(7) * t43 + t32;
t27 = t49 * t20 - t36 * t23;
t7 = t22 * t37;
t9 = t25 * t37;
t26 = g(1) * t9 + g(2) * t7 + g(3) * (t24 * pkin(7) + t40);
t5 = -t21 * t46 + t41;
t4 = t21 * t45 + t42;
t3 = t21 * t42 + t45;
t2 = -t21 * t41 + t46;
t1 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * (g(1) * (t25 * rSges(3,3) + t17) + g(2) * (rSges(3,1) * t43 - rSges(3,2) * t47 + t39) + (g(1) * (-pkin(1) - t31) + g(2) * rSges(3,3)) * t22) - m(4) * (g(1) * (t25 * rSges(4,1) + t17) + g(2) * (-rSges(4,2) * t43 + rSges(4,3) * t47 + t32) + (g(1) * (-t21 * rSges(4,3) - t16 + t33 + t44) + g(2) * rSges(4,1)) * t22) - m(5) * (g(1) * (t5 * rSges(5,1) - t4 * rSges(5,2) + t38) + g(2) * (t3 * rSges(5,1) - t2 * rSges(5,2) + rSges(5,3) * t43 + t28) + (t34 * t24 + t33) * t53) - m(6) * (g(1) * (t36 * t4 + t49 * t5 + t38) + g(2) * (rSges(6,2) * t43 + t36 * t2 + t49 * t3 + t28) + (t35 * t24 + t33) * t53), -m(3) * (g(3) * t31 + t55 * (-rSges(3,1) * t21 - rSges(3,2) * t24)) - m(4) * (g(1) * (rSges(4,3) * t43 + t9) + g(2) * (t22 * t24 * rSges(4,3) + t7) + g(3) * (t40 - t44) + (g(3) * rSges(4,3) + t55 * (rSges(4,2) - pkin(2))) * t21) - m(5) * ((g(3) * rSges(5,3) + t55 * t29) * t24 + (g(3) * t29 + t55 * t34) * t21 + t26) - m(6) * ((g(3) * rSges(6,2) + t55 * t27) * t24 + (g(3) * t27 + t55 * t35) * t21 + t26), (-m(4) - m(5) - m(6)) * (t55 * t21 - t50), -m(5) * (g(1) * (-t2 * rSges(5,1) - t3 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t5 * rSges(5,2))) - m(6) * (g(1) * (-t49 * t2 + t36 * t3) + g(2) * (-t36 * t5 + t49 * t4)) + (-m(5) * (-rSges(5,1) * t23 + rSges(5,2) * t20) - m(6) * (-t36 * t20 - t49 * t23)) * t50, -m(6) * (g(1) * t2 - g(2) * t4 + t23 * t50)];
taug = t1(:);
