% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR16_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t22 = t16 * pkin(3) - t19 * qJ(4);
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t34 = -g(1) * t17 + g(2) * t20;
t1 = g(3) * t19 - t16 * t34;
t31 = g(3) * t16;
t28 = t17 * t19;
t15 = sin(qJ(5));
t27 = t20 * t15;
t18 = cos(qJ(5));
t26 = t20 * t18;
t25 = t20 * pkin(1) + t17 * qJ(2);
t10 = g(1) * t20 + g(2) * t17;
t12 = t20 * qJ(2);
t8 = t10 * t19;
t7 = t10 * t16;
t6 = -t15 * t28 + t26;
t5 = -t18 * t28 - t27;
t4 = -t17 * t18 - t19 * t27;
t3 = t17 * t15 - t19 * t26;
t2 = -t19 * t34 - t31;
t9 = [0, -t34, t10, t34, -t10, -g(1) * (-t17 * pkin(1) + t12) - g(2) * t25, 0, 0, 0, 0, 0, -t7, -t8, -t34, t7, t8, -g(1) * (t22 * t20 + t12) - g(2) * (t20 * pkin(6) + t25) + (-g(1) * (-pkin(1) - pkin(6)) - g(2) * t22) * t17, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, t2, -t1, g(3) * t22 + t34 * (pkin(3) * t19 + qJ(4) * t16), 0, 0, 0, 0, 0, -t1 * t15, -t1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 - t18 * t31, g(1) * t6 - g(2) * t4 + t15 * t31;];
taug_reg = t9;
