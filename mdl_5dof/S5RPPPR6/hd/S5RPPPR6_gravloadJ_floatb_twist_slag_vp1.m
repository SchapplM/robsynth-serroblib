% Calculate Gravitation load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (121->78), mult. (245->108), div. (0->0), fcn. (245->8), ass. (0->35)
t46 = rSges(6,3) + pkin(6);
t45 = -m(5) - m(6);
t25 = sin(qJ(1));
t44 = g(1) * t25;
t20 = sin(pkin(8));
t43 = t20 * t25;
t21 = sin(pkin(7));
t27 = cos(qJ(1));
t42 = t21 * t27;
t22 = cos(pkin(8));
t41 = t22 * t27;
t23 = cos(pkin(7));
t24 = sin(qJ(5));
t40 = t23 * t24;
t26 = cos(qJ(5));
t39 = t23 * t26;
t38 = t23 * t27;
t37 = t25 * t22;
t36 = -pkin(2) - qJ(4);
t35 = t27 * pkin(1) + t25 * qJ(2);
t16 = t27 * qJ(2);
t34 = t27 * pkin(3) + t16;
t33 = t21 * qJ(3);
t32 = -m(4) + t45;
t31 = -pkin(1) - t33;
t30 = pkin(2) * t38 + t27 * t33 + t35;
t29 = g(1) * t27 + g(2) * t25;
t28 = t25 * pkin(3) + qJ(4) * t38 + t30;
t8 = -t21 * t43 + t41;
t7 = t20 * t27 + t21 * t37;
t6 = t20 * t42 + t37;
t5 = -t21 * t41 + t43;
t2 = t24 * t38 + t6 * t26;
t1 = -t6 * t24 + t26 * t38;
t3 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t25 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t27 + t16) + g(2) * (rSges(3,1) * t38 - rSges(3,2) * t42 + t35) + (g(1) * (-rSges(3,1) * t23 + rSges(3,2) * t21 - pkin(1)) + g(2) * rSges(3,3)) * t25) - m(4) * (g(1) * (rSges(4,1) * t27 + t16) + g(2) * (-rSges(4,2) * t38 + rSges(4,3) * t42 + t30) + (g(1) * (-rSges(4,3) * t21 + t31 + (rSges(4,2) - pkin(2)) * t23) + g(2) * rSges(4,1)) * t25) - m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t7 + t34) + g(2) * (t6 * rSges(5,1) - t5 * rSges(5,2) + rSges(5,3) * t38 + t28) + ((-rSges(5,3) + t36) * t23 + t31) * t44) - m(6) * (g(2) * (rSges(6,1) * t2 + rSges(6,2) * t1 + pkin(4) * t6 + t46 * t5 + t28) + ((-rSges(6,1) * t24 - rSges(6,2) * t26 + t36) * t23 + t31) * t44 + (t34 + (t26 * rSges(6,1) - t24 * rSges(6,2) + pkin(4)) * t8 + t46 * t7) * g(1)), (-m(3) + t32) * (-g(2) * t27 + t44), t32 * (-g(3) * t23 + t29 * t21), t45 * (g(3) * t21 + t29 * t23), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * ((t24 * t8 + t25 * t39) * rSges(6,1) + (-t25 * t40 + t26 * t8) * rSges(6,2)) + g(3) * ((t20 * t40 + t21 * t26) * rSges(6,1) + (t20 * t39 - t21 * t24) * rSges(6,2)))];
taug = t3(:);
