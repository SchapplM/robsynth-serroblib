% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t23 = t19 * pkin(4) + t17 * qJ(5);
t32 = -pkin(3) - t23;
t16 = qJ(1) + pkin(8);
t15 = qJ(3) + t16;
t14 = cos(t15);
t13 = sin(t15);
t31 = g(1) * t13;
t5 = -g(2) * t14 + t31;
t6 = g(1) * t14 + g(2) * t13;
t25 = t13 * pkin(7) - t32 * t14;
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t24 = g(1) * t18 - g(2) * t20;
t21 = t32 * t31;
t11 = t14 * pkin(7);
t4 = t5 * t19;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t19;
t1 = -g(3) * t19 + t6 * t17;
t7 = [0, t24, g(1) * t20 + g(2) * t18, t24 * pkin(1), 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-pkin(2) * sin(t16) - t18 * pkin(1) + t11) - g(2) * (pkin(2) * cos(t16) + t20 * pkin(1) + t25) - t21; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t11 - g(2) * t25 - t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t23 + t6 * (pkin(4) * t17 - qJ(5) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
