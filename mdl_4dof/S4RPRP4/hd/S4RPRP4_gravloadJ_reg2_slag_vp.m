% Calculate inertial parameters regressor of gravitation load for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t12 = qJ(1) + pkin(6);
t10 = cos(t12);
t9 = sin(t12);
t5 = g(1) * t10 + g(2) * t9;
t25 = g(1) * t9;
t16 = cos(qJ(1));
t22 = t16 * pkin(1) + t10 * pkin(2) + t9 * pkin(5);
t14 = sin(qJ(1));
t21 = -t14 * pkin(1) + t10 * pkin(5);
t20 = -g(2) * t10 + t25;
t19 = g(1) * t14 - g(2) * t16;
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t18 = t15 * pkin(3) + t13 * qJ(4);
t4 = t20 * t15;
t3 = t20 * t13;
t2 = g(3) * t13 + t5 * t15;
t1 = -g(3) * t15 + t5 * t13;
t6 = [0, 0, 0, 0, 0, 0, t19, g(1) * t16 + g(2) * t14, 0, 0, 0, 0, 0, 0, 0, 0, t20, t5, 0, t19 * pkin(1), 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-t9 * pkin(2) + t21) - g(2) * t22, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * t21 - g(2) * (t18 * t10 + t22) - (-pkin(2) - t18) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t18 + t5 * (pkin(3) * t13 - qJ(4) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t6;
