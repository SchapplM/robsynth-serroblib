% Calculate minimal parameter regressor of gravitation load for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:32:08
% EndTime: 2019-05-05 14:32:09
% DurationCPUTime: 0.20s
% Computational Cost: add. (155->53), mult. (191->70), div. (0->0), fcn. (194->10), ass. (0->36)
t25 = sin(qJ(1));
t26 = cos(qJ(1));
t45 = -g(1) * t25 + g(2) * t26;
t19 = pkin(9) + qJ(4);
t11 = sin(t19);
t13 = cos(t19);
t2 = -g(3) * t11 - t13 * t45;
t41 = g(3) * t13;
t18 = pkin(10) + qJ(6);
t10 = sin(t18);
t40 = t25 * t10;
t12 = cos(t18);
t39 = t25 * t12;
t20 = sin(pkin(10));
t38 = t25 * t20;
t22 = cos(pkin(10));
t37 = t25 * t22;
t36 = t26 * t10;
t35 = t26 * t12;
t34 = t26 * t20;
t33 = t26 * t22;
t32 = t26 * pkin(1) + t25 * qJ(2);
t31 = g(2) * t32;
t9 = g(1) * t26 + g(2) * t25;
t29 = -t11 * pkin(4) + t13 * qJ(5);
t21 = sin(pkin(9));
t28 = pkin(3) * t21 - t29;
t24 = -pkin(7) - qJ(3);
t15 = t26 * qJ(2);
t7 = t9 * t13;
t6 = t11 * t35 - t40;
t5 = t11 * t36 + t39;
t4 = t11 * t39 + t36;
t3 = -t11 * t40 + t35;
t1 = -t11 * t45 + t41;
t8 = [0, -t45, t9, t45, -t9, -g(1) * (-t25 * pkin(1) + t15) - t31, -t9 * t21, -t9 * cos(pkin(9)) -t45, -g(1) * (t15 + (-pkin(1) - qJ(3)) * t25) - g(2) * (t26 * qJ(3) + t32) 0, 0, 0, 0, 0, -t9 * t11, -t7, -g(1) * (t11 * t33 - t38) - g(2) * (t11 * t37 + t34) -g(1) * (-t11 * t34 - t37) - g(2) * (-t11 * t38 + t33) t7, -g(1) * t15 - t31 + (-g(1) * t28 + g(2) * t24) * t26 + (-g(1) * (-pkin(1) + t24) - g(2) * t28) * t25, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, t45, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2 * t22, t2 * t20, -t1, -g(3) * t29 + t45 * (pkin(4) * t13 + qJ(5) * t11) 0, 0, 0, 0, 0, -t2 * t12, t2 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t10 * t41, g(1) * t4 - g(2) * t6 + t12 * t41;];
taug_reg  = t8;
