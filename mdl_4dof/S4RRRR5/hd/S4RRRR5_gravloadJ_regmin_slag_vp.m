% Calculate minimal parameter regressor of gravitation load for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = sin(qJ(2));
t18 = cos(qJ(2));
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t23 = g(1) * t19 + g(2) * t16;
t21 = -g(3) * t18 + t23 * t15;
t30 = g(3) * t15;
t28 = t16 * t18;
t13 = qJ(3) + qJ(4);
t11 = sin(t13);
t27 = t19 * t11;
t12 = cos(t13);
t26 = t19 * t12;
t14 = sin(qJ(3));
t25 = t19 * t14;
t17 = cos(qJ(3));
t24 = t19 * t17;
t22 = g(1) * t16 - g(2) * t19;
t10 = t16 * t14 + t18 * t24;
t9 = t16 * t17 - t18 * t25;
t8 = -t17 * t28 + t25;
t7 = t14 * t28 + t24;
t6 = t16 * t11 + t18 * t26;
t5 = t16 * t12 - t18 * t27;
t4 = -t12 * t28 + t27;
t3 = t11 * t28 + t26;
t2 = g(1) * t6 - g(2) * t4 + t12 * t30;
t1 = -g(1) * t5 + g(2) * t3 + t11 * t30;
t20 = [0, t22, t23, 0, 0, 0, 0, 0, t22 * t18, -t22 * t15, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t21, t23 * t18 + t30, 0, 0, 0, 0, 0, t21 * t17, -t21 * t14, 0, 0, 0, 0, 0, t21 * t12, -t21 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t14 * t30, g(1) * t10 - g(2) * t8 + t17 * t30, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t20;
