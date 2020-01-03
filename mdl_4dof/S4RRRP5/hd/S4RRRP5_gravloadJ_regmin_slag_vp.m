% Calculate minimal parameter regressor of gravitation load for
% S4RRRP5
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
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = cos(qJ(2));
t14 = qJ(2) + qJ(3);
t11 = sin(t14);
t12 = cos(t14);
t25 = t12 * pkin(3) + t11 * qJ(4);
t26 = t17 * pkin(2) + t25;
t24 = pkin(3) * t11;
t23 = qJ(4) * t12;
t15 = sin(qJ(2));
t22 = -pkin(2) * t15 - t24;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t5 = g(1) * t18 + g(2) * t16;
t21 = g(1) * t16 - g(2) * t18;
t20 = pkin(1) + t26;
t19 = -pkin(6) - pkin(5);
t7 = t18 * t23;
t6 = t16 * t23;
t4 = t21 * t12;
t3 = t21 * t11;
t2 = g(3) * t11 + t5 * t12;
t1 = -g(3) * t12 + t5 * t11;
t8 = [0, t21, t5, 0, 0, 0, 0, 0, t21 * t17, -t21 * t15, 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3, (g(1) * t19 - g(2) * t20) * t18 + (g(1) * t20 + g(2) * t19) * t16; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t17 + t5 * t15, g(3) * t15 + t5 * t17, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (t22 * t18 + t7) - g(2) * (t22 * t16 + t6) - g(3) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t18 * t24 + t7) - g(2) * (-t16 * t24 + t6) - g(3) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t8;
