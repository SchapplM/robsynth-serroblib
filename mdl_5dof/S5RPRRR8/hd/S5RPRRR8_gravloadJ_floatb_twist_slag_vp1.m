% Calculate Gravitation load on the joints for
% S5RPRRR8
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:42
% EndTime: 2019-12-31 19:05:43
% DurationCPUTime: 0.33s
% Computational Cost: add. (193->56), mult. (305->76), div. (0->0), fcn. (327->8), ass. (0->30)
t17 = cos(qJ(4));
t15 = qJ(4) + qJ(5);
t10 = cos(t15);
t9 = sin(t15);
t29 = -t10 * rSges(6,1) + t9 * rSges(6,2);
t49 = -t17 * pkin(4) + t29;
t34 = sin(qJ(3));
t35 = sin(qJ(1));
t36 = cos(qJ(3));
t37 = cos(qJ(1));
t1 = -t35 * t34 - t37 * t36;
t2 = t37 * t34 - t35 * t36;
t23 = pkin(3) - t49;
t32 = rSges(6,3) + pkin(8) + pkin(7);
t48 = -(g(1) * t32 - g(2) * t23) * t1 - (g(1) * t23 + g(2) * t32) * t2;
t16 = sin(qJ(4));
t24 = -t17 * rSges(5,1) + t16 * rSges(5,2);
t22 = pkin(3) - t24;
t38 = rSges(5,3) + pkin(7);
t47 = -(g(1) * t38 - g(2) * t22) * t1 - (g(1) * t22 + g(2) * t38) * t2;
t46 = g(1) * t1 + g(2) * t2;
t31 = t37 * pkin(1) + t35 * qJ(2);
t30 = t37 * pkin(2) + t31;
t28 = -t35 * pkin(1) + t37 * qJ(2);
t27 = -t2 * rSges(4,1) + t1 * rSges(4,2);
t26 = t1 * rSges(4,1) + t2 * rSges(4,2);
t25 = rSges(6,1) * t9 + rSges(6,2) * t10;
t21 = -t35 * pkin(2) + t28;
t20 = g(1) * t21 + g(2) * t30;
t3 = [-m(2) * (g(1) * (-t35 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) - t35 * rSges(2,2))) - m(3) * (g(1) * (-t35 * rSges(3,1) + t37 * rSges(3,3) + t28) + g(2) * (t37 * rSges(3,1) + t35 * rSges(3,3) + t31)) - m(4) * (g(1) * (t21 - t27) + g(2) * (-t26 + t30)) - m(5) * (t20 - t47) - m(6) * (t20 - t48), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t35 - g(2) * t37), -m(4) * (g(1) * t27 + g(2) * t26) - m(5) * t47 - m(6) * t48, (-m(5) * t24 - m(6) * t49) * g(3) + t46 * (-m(5) * (rSges(5,1) * t16 + rSges(5,2) * t17) - m(6) * (pkin(4) * t16 + t25)), -m(6) * (g(3) * t29 + t46 * t25)];
taug = t3(:);
