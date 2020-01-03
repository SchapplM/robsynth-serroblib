% Calculate inertial parameters regressor of gravitation load for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = sin(qJ(1));
t23 = cos(qJ(1));
t24 = t23 * pkin(1) + t22 * qJ(2);
t21 = cos(pkin(6));
t20 = sin(pkin(6));
t19 = t23 * pkin(2) + t24;
t1 = -t22 * t20 - t23 * t21;
t2 = t23 * t20 - t22 * t21;
t18 = g(1) * t2 - g(2) * t1;
t17 = g(1) * t1 + g(2) * t2;
t16 = -t22 * pkin(1) + t23 * qJ(2);
t15 = -t22 * pkin(2) + t16;
t14 = cos(qJ(4));
t13 = sin(qJ(4));
t4 = g(1) * t23 + g(2) * t22;
t3 = g(1) * t22 - g(2) * t23;
t5 = [0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, -t4, -g(1) * t16 - g(2) * t24, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -g(1) * t15 - g(2) * t19, 0, 0, 0, 0, 0, 0, -t18 * t14, t18 * t13, -t17, -g(1) * (t2 * pkin(3) + t1 * pkin(5) + t15) - g(2) * (-t1 * pkin(3) + t2 * pkin(5) + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t14 - t17 * t13, -g(3) * t13 - t17 * t14, 0, 0;];
taug_reg = t5;
