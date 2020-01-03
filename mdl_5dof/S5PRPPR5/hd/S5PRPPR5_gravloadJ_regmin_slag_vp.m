% Calculate minimal parameter regressor of gravitation load for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = sin(pkin(7));
t17 = cos(pkin(7));
t24 = g(1) * t17 + g(2) * t15;
t19 = sin(qJ(2));
t32 = t24 * t19;
t14 = sin(pkin(8));
t16 = cos(pkin(8));
t21 = cos(qJ(2));
t23 = t21 * t14 - t19 * t16;
t31 = g(3) * t23;
t29 = pkin(2) * t19;
t26 = t21 * pkin(2) + t19 * qJ(3);
t25 = qJ(3) * t21;
t8 = t19 * t14 + t21 * t16;
t22 = g(3) * t8 + t24 * t23;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t10 = t17 * t25;
t9 = t15 * t25;
t6 = t8 * t17;
t4 = t8 * t15;
t2 = g(3) * t19 + t24 * t21;
t1 = -g(3) * t21 + t32;
t3 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t1, t2, t1, -t2, -g(1) * (-t17 * t29 + t10) - g(2) * (-t15 * t29 + t9) - g(3) * t26, -t22, -g(1) * t6 - g(2) * t4 + t31, -g(1) * t10 - g(2) * t9 - g(3) * (t21 * pkin(3) + t26) + (pkin(2) + pkin(3)) * t32, 0, 0, 0, 0, 0, -t22 * t20, t22 * t18; 0, 0, 0, 0, 0, 0, -t1, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t20 - t6 * t18) - g(2) * (t17 * t20 - t4 * t18) - t18 * t31, -g(1) * (t15 * t18 - t6 * t20) - g(2) * (-t17 * t18 - t4 * t20) - t20 * t31;];
taug_reg = t3;
