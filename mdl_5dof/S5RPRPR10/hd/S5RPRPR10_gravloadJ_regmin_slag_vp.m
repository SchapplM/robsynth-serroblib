% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(qJ(5));
t28 = qJ(3) + pkin(8);
t12 = sin(t28);
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t27 = cos(t28);
t1 = t22 * t12 - t19 * t27;
t2 = t19 * t12 + t22 * t27;
t25 = g(1) * t1 + g(2) * t2;
t33 = t25 * t17;
t20 = cos(qJ(5));
t32 = t25 * t20;
t18 = sin(qJ(3));
t31 = t19 * t18;
t30 = t22 * t18;
t29 = t22 * pkin(1) + t19 * qJ(2);
t26 = g(1) * t2 - g(2) * t1;
t21 = cos(qJ(3));
t3 = -t22 * t21 - t31;
t4 = -t19 * t21 + t30;
t24 = g(1) * t4 - g(2) * t3;
t23 = g(1) * t3 + g(2) * t4;
t14 = t22 * qJ(2);
t11 = t21 * pkin(3) + pkin(2);
t6 = g(1) * t22 + g(2) * t19;
t5 = g(1) * t19 - g(2) * t22;
t7 = [0, t5, t6, t5, -t6, -g(1) * (-t19 * pkin(1) + t14) - g(2) * t29, 0, -t24, t23, -g(1) * (pkin(3) * t30 + t14 + (-pkin(1) - t11) * t19) - g(2) * (pkin(3) * t31 + t22 * t11 + t29), 0, 0, 0, 0, 0, -t32, t33; 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t24, -t23, t24 * pkin(3), 0, 0, 0, 0, 0, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t20 + t26 * t17, -g(3) * t17 + t26 * t20;];
taug_reg = t7;
