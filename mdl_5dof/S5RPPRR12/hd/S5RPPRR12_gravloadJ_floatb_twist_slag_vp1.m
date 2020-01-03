% Calculate Gravitation load on the joints for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:06:58
% DurationCPUTime: 0.35s
% Computational Cost: add. (146->68), mult. (195->94), div. (0->0), fcn. (169->8), ass. (0->33)
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t50 = -g(1) * t20 + g(2) * t22;
t47 = g(1) * t22 + g(2) * t20;
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t23 = m(5) * rSges(5,1) + m(6) * (rSges(6,1) * t21 - rSges(6,2) * t19 + pkin(4));
t37 = rSges(6,3) + pkin(7);
t46 = -m(5) * rSges(5,2) + m(6) * t37;
t16 = sin(pkin(8));
t43 = pkin(3) * t16;
t15 = pkin(8) + qJ(4);
t10 = sin(t15);
t42 = pkin(4) * t10;
t36 = t19 * t22;
t35 = t20 * t19;
t34 = t20 * t21;
t33 = t21 * t22;
t32 = t22 * pkin(1) + t20 * qJ(2);
t31 = rSges(4,3) + qJ(3);
t30 = -m(4) - m(5) - m(6);
t13 = t22 * qJ(2);
t18 = -pkin(6) - qJ(3);
t29 = t20 * t18 + t22 * t43 + t13;
t28 = t20 * t43 + t32;
t27 = rSges(4,1) * t16 + rSges(4,2) * cos(pkin(8));
t11 = cos(t15);
t26 = rSges(5,1) * t10 + rSges(5,2) * t11;
t4 = t10 * t33 - t35;
t3 = t10 * t36 + t34;
t2 = t10 * t34 + t36;
t1 = -t10 * t35 + t33;
t5 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - rSges(2,2) * t22) + g(2) * (rSges(2,1) * t22 - t20 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t22 + t13 + (rSges(3,2) - pkin(1)) * t20) + g(2) * (-rSges(3,2) * t22 + t20 * rSges(3,3) + t32)) - m(4) * (g(1) * t13 + g(2) * t32 + (g(1) * t27 + g(2) * t31) * t22 + (g(1) * (-pkin(1) - t31) + g(2) * t27) * t20) - m(5) * (g(1) * t29 + g(2) * t28 + (g(1) * t26 + g(2) * (rSges(5,3) - t18)) * t22 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t26) * t20) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t20 * pkin(1) + t22 * t42 + t29) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t18 * t22 + t20 * t42 + t28) - t47 * t11 * t37), -(-m(3) + t30) * t50, t30 * t47, (t23 * t10 - t11 * t46) * g(3) + t50 * (t46 * t10 + t23 * t11), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t19 - rSges(6,2) * t21) * t11)];
taug = t5(:);
