% Calculate minimal parameter regressor of gravitation load for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(qJ(3));
t19 = qJ(1) + pkin(7);
t13 = sin(t19);
t14 = cos(t19);
t5 = g(1) * t14 + g(2) * t13;
t38 = t5 * t20;
t15 = t20 * qJ(4);
t22 = cos(qJ(3));
t30 = t22 * pkin(3) + t15;
t36 = pkin(3) * t20;
t35 = g(1) * t13;
t32 = t22 * pkin(4);
t31 = t14 * t22;
t29 = qJ(4) * t22;
t21 = sin(qJ(1));
t28 = -t21 * pkin(1) + t14 * pkin(6);
t23 = cos(qJ(1));
t27 = t23 * pkin(1) + pkin(3) * t31 + t13 * pkin(6) + (pkin(2) + t15) * t14;
t26 = -g(2) * t14 + t35;
t25 = g(1) * t21 - g(2) * t23;
t24 = -pkin(2) - t30;
t8 = t14 * t29;
t6 = t13 * t29;
t4 = t26 * t22;
t3 = t26 * t20;
t2 = g(3) * t20 + t5 * t22;
t1 = -g(3) * t22 + t38;
t7 = [0, t25, g(1) * t23 + g(2) * t21, t25 * pkin(1), 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3, -g(1) * t28 - g(2) * t27 - t24 * t35, t4, t3, t5, -g(1) * (-t14 * qJ(5) + t28) - g(2) * (pkin(4) * t31 + t27) + (-g(1) * (t24 - t32) + g(2) * qJ(5)) * t13; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t14 * t36 + t8) - g(2) * (-t13 * t36 + t6) - g(3) * t30, t1, -t2, 0, -g(1) * t8 - g(2) * t6 - g(3) * (t30 + t32) + (pkin(3) + pkin(4)) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26;];
taug_reg = t7;
