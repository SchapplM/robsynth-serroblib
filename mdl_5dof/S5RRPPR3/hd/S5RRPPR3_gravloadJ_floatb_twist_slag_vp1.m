% Calculate Gravitation load on the joints for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (221->56), mult. (138->65), div. (0->0), fcn. (102->8), ass. (0->31)
t49 = rSges(6,3) + pkin(7);
t24 = sin(qJ(5));
t26 = cos(qJ(5));
t48 = rSges(6,1) * t24 + rSges(6,2) * t26;
t23 = qJ(1) + qJ(2);
t19 = pkin(8) + t23;
t17 = cos(t19);
t16 = sin(t19);
t43 = g(1) * t16;
t47 = -g(2) * t17 + t43;
t46 = -m(5) - m(6);
t25 = sin(qJ(1));
t45 = pkin(1) * t25;
t20 = sin(t23);
t44 = pkin(2) * t20;
t21 = cos(t23);
t18 = pkin(2) * t21;
t39 = t17 * pkin(3) + t16 * qJ(4) + t18;
t38 = t17 * qJ(4) - t44;
t37 = t21 * rSges(3,1) - rSges(3,2) * t20;
t36 = t17 * rSges(4,1) - rSges(4,2) * t16 + t18;
t35 = t48 * t17 + t38;
t34 = (-pkin(3) - t49) * t43;
t33 = -rSges(3,1) * t20 - rSges(3,2) * t21;
t31 = -rSges(5,2) * t17 + t16 * rSges(5,3) + t39;
t30 = t48 * t16 + t49 * t17 + t39;
t29 = t17 * rSges(5,3) + t38 + (rSges(5,2) - pkin(3)) * t16;
t28 = -rSges(4,1) * t16 - rSges(4,2) * t17 - t44;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t25 * rSges(2,2))) - m(3) * (g(1) * (t33 - t45) + g(2) * (t22 + t37)) - m(4) * (g(1) * (t28 - t45) + g(2) * (t22 + t36)) - m(5) * (g(1) * (t29 - t45) + g(2) * (t22 + t31)) - m(6) * (g(1) * (t35 - t45) + g(2) * (t22 + t30) + t34), -m(3) * (g(1) * t33 + g(2) * t37) - m(4) * (g(1) * t28 + g(2) * t36) - m(5) * (g(1) * t29 + g(2) * t31) - m(6) * (g(1) * t35 + g(2) * t30 + t34), (-m(4) + t46) * g(3), t46 * t47, -m(6) * (-g(3) * t48 + t47 * (rSges(6,1) * t26 - rSges(6,2) * t24))];
taug = t1(:);
