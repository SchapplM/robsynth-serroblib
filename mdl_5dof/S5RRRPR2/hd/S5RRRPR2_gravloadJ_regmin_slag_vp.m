% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:52
% EndTime: 2022-01-20 11:30:53
% DurationCPUTime: 0.09s
% Computational Cost: add. (128->23), mult. (76->33), div. (0->0), fcn. (70->10), ass. (0->24)
t17 = qJ(1) + qJ(2);
t16 = qJ(3) + t17;
t13 = cos(t16);
t15 = cos(t17);
t25 = pkin(2) * t15 + pkin(3) * t13;
t11 = pkin(9) + t16;
t7 = sin(t11);
t8 = cos(t11);
t24 = g(1) * t8 + g(2) * t7;
t23 = g(1) * t7 - g(2) * t8;
t12 = sin(t16);
t14 = sin(t17);
t22 = -pkin(2) * t14 - pkin(3) * t12;
t3 = g(1) * t12 - g(2) * t13;
t21 = cos(qJ(1));
t20 = cos(qJ(5));
t19 = sin(qJ(1));
t18 = sin(qJ(5));
t6 = g(1) * t15 + g(2) * t14;
t5 = g(1) * t14 - g(2) * t15;
t4 = g(1) * t13 + g(2) * t12;
t2 = t23 * t20;
t1 = t23 * t18;
t9 = [0, g(1) * t19 - g(2) * t21, g(1) * t21 + g(2) * t19, 0, t5, t6, 0, t3, t4, -g(1) * (-t19 * pkin(1) + t22) - g(2) * (t21 * pkin(1) + t25), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t5, t6, 0, t3, t4, -g(1) * t22 - g(2) * t25, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, t3, t4, t3 * pkin(3), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t20 + t24 * t18, g(3) * t18 + t24 * t20;];
taug_reg = t9;
