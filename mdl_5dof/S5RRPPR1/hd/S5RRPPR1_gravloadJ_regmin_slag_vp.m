% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:42
% EndTime: 2022-01-20 09:51:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (127->29), mult. (80->34), div. (0->0), fcn. (72->9), ass. (0->25)
t20 = qJ(1) + qJ(2);
t16 = sin(t20);
t29 = pkin(2) * t16;
t22 = sin(qJ(1));
t28 = t22 * pkin(1);
t15 = pkin(8) + t20;
t10 = sin(t15);
t11 = cos(t15);
t17 = cos(t20);
t12 = pkin(2) * t17;
t27 = t11 * pkin(3) + t10 * qJ(4) + t12;
t26 = g(1) * t11 + g(2) * t10;
t25 = g(1) * t10 - g(2) * t11;
t5 = g(1) * t16 - g(2) * t17;
t24 = -t10 * pkin(3) + t11 * qJ(4) - t29;
t23 = cos(qJ(1));
t19 = pkin(9) + qJ(5);
t18 = t23 * pkin(1);
t14 = cos(t19);
t13 = sin(t19);
t6 = g(1) * t17 + g(2) * t16;
t3 = t25 * cos(pkin(9));
t2 = t25 * t14;
t1 = t25 * t13;
t4 = [0, g(1) * t22 - g(2) * t23, g(1) * t23 + g(2) * t22, 0, t5, t6, -g(1) * (-t28 - t29) - g(2) * (t12 + t18), t3, -t26, -g(1) * (t24 - t28) - g(2) * (t18 + t27), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t5, t6, t5 * pkin(2), t3, -t26, -g(1) * t24 - g(2) * t27, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t14 + t26 * t13, g(3) * t13 + t26 * t14;];
taug_reg = t4;
