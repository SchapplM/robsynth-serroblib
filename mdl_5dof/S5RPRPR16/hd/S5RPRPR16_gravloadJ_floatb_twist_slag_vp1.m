% Calculate Gravitation load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR16_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:46
% DurationCPUTime: 0.34s
% Computational Cost: add. (120->86), mult. (233->118), div. (0->0), fcn. (204->6), ass. (0->36)
t47 = rSges(6,3) + pkin(7);
t20 = cos(qJ(3));
t12 = t20 * qJ(4);
t46 = -rSges(5,3) * t20 - t12;
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t24 = rSges(6,1) * t16 + rSges(6,2) * t19;
t45 = -m(5) - m(6);
t44 = -pkin(1) - pkin(6);
t18 = sin(qJ(1));
t43 = g(1) * t18;
t21 = cos(qJ(1));
t42 = g(2) * t21;
t17 = sin(qJ(3));
t41 = g(3) * t17;
t40 = rSges(5,2) - pkin(3);
t36 = t17 * t18;
t35 = t17 * t21;
t34 = t18 * t20;
t33 = t20 * t21;
t13 = t21 * qJ(2);
t32 = pkin(3) * t35 + t13;
t31 = t21 * pkin(1) + t18 * qJ(2);
t30 = -pkin(3) - t47;
t29 = t21 * pkin(6) + t31;
t28 = g(1) * (pkin(3) * t34 + qJ(4) * t36);
t27 = pkin(3) * t36 + t29;
t26 = -t42 + t43;
t25 = rSges(4,1) * t17 + rSges(4,2) * t20;
t23 = t47 * t17 - t12;
t22 = -rSges(5,2) * t17 + t46;
t5 = -t16 * t34 + t19 * t21;
t4 = -t16 * t21 - t19 * t34;
t3 = -t16 * t33 - t18 * t19;
t2 = t18 * t16 - t19 * t33;
t1 = [-m(2) * (g(1) * (-t18 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t18 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t21 + t13 + (rSges(3,2) - pkin(1)) * t18) + g(2) * (-rSges(3,2) * t21 + t18 * rSges(3,3) + t31)) - m(4) * (g(1) * (rSges(4,1) * t35 + rSges(4,2) * t33 + t13) + g(2) * (rSges(4,3) * t21 + t29) + (g(1) * (-rSges(4,3) + t44) + g(2) * t25) * t18) - m(5) * (g(1) * t32 + g(2) * t27 + (g(2) * rSges(5,1) + g(1) * t22) * t21 + (g(1) * (-rSges(5,1) + t44) + g(2) * t22) * t18) - m(6) * (g(1) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t32) + g(2) * (t5 * rSges(6,1) + t4 * rSges(6,2) + t27) + (g(2) * pkin(4) + g(1) * t23) * t21 + (g(1) * (-pkin(4) + t44) + g(2) * t23) * t18), (-m(3) - m(4) + t45) * t26, m(4) * g(3) * t25 - m(5) * (t28 + g(3) * (t17 * t40 - t46)) - m(6) * (t28 + g(3) * (t17 * t30 + t20 * t24 + t12)) + (-m(4) * (rSges(4,1) * t20 - rSges(4,2) * t17) - m(5) * (-rSges(5,2) * t20 + rSges(5,3) * t17) - m(6) * (t24 * t17 + t47 * t20)) * t43 + ((m(4) * rSges(4,1) - m(5) * t40 - m(6) * t30) * t20 + (-m(4) * rSges(4,2) - m(5) * (-rSges(5,3) - qJ(4)) - m(6) * (-qJ(4) - t24)) * t17) * t42, t45 * (-t20 * t26 + t41), -m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 + rSges(6,2) * t3) + (rSges(6,1) * t19 - rSges(6,2) * t16) * t41)];
taug = t1(:);
