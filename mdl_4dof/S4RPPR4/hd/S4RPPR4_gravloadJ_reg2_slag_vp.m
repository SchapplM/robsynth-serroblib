% Calculate inertial parameters regressor of gravitation load for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t14 = cos(qJ(1));
t10 = qJ(1) + pkin(6);
t7 = sin(t10);
t8 = cos(t10);
t17 = t14 * pkin(1) + t8 * pkin(2) + t7 * qJ(3);
t12 = sin(qJ(1));
t16 = -t12 * pkin(1) + t8 * qJ(3);
t2 = g(1) * t8 + g(2) * t7;
t1 = g(1) * t7 - g(2) * t8;
t15 = g(1) * t12 - g(2) * t14;
t13 = cos(qJ(4));
t11 = sin(qJ(4));
t3 = [0, 0, 0, 0, 0, 0, t15, g(1) * t14 + g(2) * t12, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t15 * pkin(1), 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t7 * pkin(2) + t16) - g(2) * t17, 0, 0, 0, 0, 0, 0, -t2 * t11, -t2 * t13, t1, -g(1) * ((-pkin(2) - pkin(5)) * t7 + t16) - g(2) * (t8 * pkin(5) + t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t11 - t1 * t13, g(3) * t13 + t1 * t11, 0, 0;];
taug_reg = t3;
