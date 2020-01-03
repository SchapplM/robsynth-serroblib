% Calculate Gravitation load on the joints for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR15_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:36:40
% DurationCPUTime: 0.44s
% Computational Cost: add. (156->74), mult. (252->102), div. (0->0), fcn. (227->8), ass. (0->33)
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t49 = -g(1) * t19 + g(2) * t21;
t34 = rSges(6,3) + pkin(7) + qJ(4);
t18 = sin(qJ(3));
t20 = cos(qJ(3));
t15 = sin(pkin(8));
t16 = cos(pkin(8));
t26 = rSges(5,1) * t16 - rSges(5,2) * t15 + pkin(3);
t32 = rSges(5,3) + qJ(4);
t46 = t26 * t18 - t32 * t20;
t7 = t16 * pkin(4) + pkin(3);
t14 = pkin(8) + qJ(5);
t8 = sin(t14);
t9 = cos(t14);
t22 = m(5) * t26 + m(6) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t7) + m(4) * rSges(4,1);
t45 = -m(4) * rSges(4,2) + m(5) * t32 + m(6) * t34;
t44 = -m(5) - m(6);
t43 = -pkin(1) - pkin(6);
t40 = pkin(4) * t15;
t37 = t19 * t18;
t36 = t20 * rSges(4,2);
t35 = t21 * t18;
t33 = t21 * pkin(1) + t19 * qJ(2);
t31 = t21 * pkin(6) + t33;
t28 = t15 * rSges(5,1) + t16 * rSges(5,2);
t24 = t18 * t7 - t34 * t20;
t11 = t21 * qJ(2);
t4 = -t19 * t8 + t9 * t35;
t3 = t19 * t9 + t8 * t35;
t2 = t21 * t8 + t9 * t37;
t1 = t21 * t9 - t8 * t37;
t5 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - t21 * rSges(2,2)) + g(2) * (t21 * rSges(2,1) - t19 * rSges(2,2))) - m(3) * (g(1) * (t21 * rSges(3,3) + t11 + (rSges(3,2) - pkin(1)) * t19) + g(2) * (-t21 * rSges(3,2) + t19 * rSges(3,3) + t33)) - m(4) * (g(1) * (rSges(4,1) * t35 + t21 * t36 + t11) + g(2) * (t21 * rSges(4,3) + t31) + (g(1) * (-rSges(4,3) + t43) + g(2) * (t18 * rSges(4,1) + t36)) * t19) - m(5) * (g(1) * t11 + g(2) * t31 + (g(1) * t46 + g(2) * t28) * t21 + (g(1) * (-t28 + t43) + t46 * g(2)) * t19) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t11) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t31) + (g(1) * t24 + g(2) * t40) * t21 + (g(1) * (-t40 + t43) + g(2) * t24) * t19), -(-m(3) - m(4) + t44) * t49, (t22 * t18 - t20 * t45) * g(3) + t49 * (t18 * t45 + t22 * t20), t44 * (g(3) * t18 + t49 * t20), -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t4 * rSges(6,2)) + g(3) * (-rSges(6,1) * t8 - rSges(6,2) * t9) * t20)];
taug = t5(:);
