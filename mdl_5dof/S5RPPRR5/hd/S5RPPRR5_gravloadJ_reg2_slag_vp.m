% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(5));
t30 = qJ(1) + pkin(8);
t26 = sin(t30);
t27 = cos(t30);
t31 = sin(qJ(4));
t32 = cos(qJ(4));
t3 = -t26 * t31 - t27 * t32;
t4 = -t26 * t32 + t27 * t31;
t23 = g(1) * t4 - g(2) * t3;
t34 = t23 * t16;
t18 = cos(qJ(5));
t33 = t23 * t18;
t19 = cos(qJ(1));
t29 = t19 * pkin(1) + t27 * pkin(2) + t26 * qJ(3);
t28 = t27 * pkin(3) + t29;
t25 = -t4 * pkin(4) - t3 * pkin(7);
t24 = t3 * pkin(4) - t4 * pkin(7);
t2 = g(1) * t3 + g(2) * t4;
t17 = sin(qJ(1));
t22 = g(1) * t17 - g(2) * t19;
t21 = -t17 * pkin(1) - t26 * pkin(2) + t27 * qJ(3);
t20 = -t26 * pkin(3) + t21;
t6 = g(1) * t27 + g(2) * t26;
t5 = g(1) * t26 - g(2) * t27;
t1 = [0, 0, 0, 0, 0, 0, t22, g(1) * t19 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t22 * pkin(1), 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * t21 - g(2) * t29, 0, 0, 0, 0, 0, 0, -t23, t2, 0, -g(1) * t20 - g(2) * t28, 0, 0, 0, 0, 0, 0, -t33, t34, -t2, -g(1) * (t20 - t25) - g(2) * (-t24 + t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t2, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t34, t2, -g(1) * t25 - g(2) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t18 - t2 * t16, -g(3) * t16 - t2 * t18, 0, 0;];
taug_reg = t1;
