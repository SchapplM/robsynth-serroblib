% Calculate inertial parameters regressor of gravitation load for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t15 = sin(qJ(4));
t24 = sin(qJ(3));
t25 = sin(qJ(1));
t26 = cos(qJ(3));
t27 = cos(qJ(1));
t3 = -t25 * t24 - t27 * t26;
t4 = t27 * t24 - t25 * t26;
t19 = g(1) * t4 - g(2) * t3;
t29 = t19 * t15;
t16 = cos(qJ(4));
t28 = t19 * t16;
t23 = t27 * pkin(1) + t25 * qJ(2);
t22 = t27 * pkin(2) + t23;
t21 = -t4 * pkin(3) - t3 * pkin(6);
t20 = t3 * pkin(3) - t4 * pkin(6);
t2 = g(1) * t3 + g(2) * t4;
t18 = -t25 * pkin(1) + t27 * qJ(2);
t17 = -t25 * pkin(2) + t18;
t6 = g(1) * t27 + g(2) * t25;
t5 = g(1) * t25 - g(2) * t27;
t1 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * t18 - g(2) * t23, 0, 0, 0, 0, 0, 0, -t19, t2, 0, -g(1) * t17 - g(2) * t22, 0, 0, 0, 0, 0, 0, -t28, t29, -t2, -g(1) * (t17 - t21) - g(2) * (-t20 + t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t2, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, t2, -g(1) * t21 - g(2) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t16 - t2 * t15, -g(3) * t15 - t2 * t16, 0, 0;];
taug_reg = t1;
