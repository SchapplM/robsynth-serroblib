% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t18 = g(1) * t13 + g(2) * t12;
t11 = qJ(2) + qJ(3);
t8 = pkin(9) + t11;
t6 = sin(t8);
t7 = cos(t8);
t25 = -g(3) * t7 + t18 * t6;
t24 = g(3) * t6;
t14 = sin(qJ(5));
t22 = t12 * t14;
t16 = cos(qJ(5));
t21 = t12 * t16;
t20 = t13 * t14;
t19 = t13 * t16;
t10 = cos(t11);
t9 = sin(t11);
t3 = -g(3) * t10 + t18 * t9;
t17 = cos(qJ(2));
t15 = sin(qJ(2));
t4 = g(3) * t9 + t18 * t10;
t2 = t25 * t16;
t1 = t25 * t14;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(3) * t17 + t18 * t15, g(3) * t15 + t18 * t17, 0, t3, t4, -g(3) * (t17 * pkin(2) + pkin(3) * t10) - t18 * (-t15 * pkin(2) - pkin(3) * t9), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, t3, t4, t3 * pkin(3), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t7 * t20 + t21) - g(2) * (-t7 * t22 - t19) + t14 * t24, -g(1) * (-t7 * t19 - t22) - g(2) * (-t7 * t21 + t20) + t16 * t24;];
taug_reg = t5;
