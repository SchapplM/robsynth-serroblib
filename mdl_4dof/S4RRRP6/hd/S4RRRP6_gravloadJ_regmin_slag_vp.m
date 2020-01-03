% Calculate minimal parameter regressor of gravitation load for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = sin(qJ(1));
t15 = cos(qJ(1));
t20 = g(1) * t15 + g(2) * t12;
t10 = sin(qJ(3));
t11 = sin(qJ(2));
t27 = g(3) * t11;
t13 = cos(qJ(3));
t23 = t15 * t13;
t14 = cos(qJ(2));
t25 = t12 * t14;
t3 = t10 * t25 + t23;
t24 = t15 * t10;
t5 = t12 * t13 - t14 * t24;
t32 = -g(1) * t5 + g(2) * t3 + t10 * t27;
t1 = -g(3) * t14 + t20 * t11;
t21 = pkin(3) * t10 + pkin(5);
t19 = g(1) * t12 - g(2) * t15;
t8 = t13 * pkin(3) + pkin(2);
t9 = -qJ(4) - pkin(6);
t18 = -t11 * t9 + t14 * t8;
t16 = -pkin(1) - t18;
t7 = t19 * t11;
t6 = t12 * t10 + t14 * t23;
t4 = -t13 * t25 + t24;
t2 = t20 * t14 + t27;
t17 = [0, t19, t20, 0, 0, 0, 0, 0, t19 * t14, -t7, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t7, (-g(1) * t21 + g(2) * t16) * t15 + (-g(1) * t16 - g(2) * t21) * t12; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t13, -t1 * t10, -t2, -g(3) * t18 + t20 * (t11 * t8 + t14 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(1) * t6 - g(2) * t4 + t13 * t27, 0, t32 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t17;
