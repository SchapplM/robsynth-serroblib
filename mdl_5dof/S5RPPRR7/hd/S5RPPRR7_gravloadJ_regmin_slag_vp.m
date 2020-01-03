% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR7
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
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t10 = sin(qJ(4));
t13 = cos(qJ(4));
t8 = qJ(1) + pkin(8);
t6 = sin(t8);
t7 = cos(t8);
t17 = g(1) * t6 - g(2) * t7;
t23 = -g(3) * t10 + t17 * t13;
t21 = g(3) * t13;
t9 = sin(qJ(5));
t20 = t10 * t9;
t12 = cos(qJ(5));
t19 = t10 * t12;
t18 = -g(1) * t7 - g(2) * t6;
t11 = sin(qJ(1));
t14 = cos(qJ(1));
t16 = g(1) * t11 - g(2) * t14;
t4 = t7 * t19 - t6 * t9;
t3 = t6 * t12 + t7 * t20;
t2 = t6 * t19 + t7 * t9;
t1 = t7 * t12 - t6 * t20;
t5 = [0, t16, g(1) * t14 + g(2) * t11, t16 * pkin(1), -t17, t18, -g(1) * (-t11 * pkin(1) - t6 * pkin(2) + t7 * qJ(3)) - g(2) * (t14 * pkin(1) + t7 * pkin(2) + t6 * qJ(3)), 0, 0, 0, 0, 0, t18 * t10, t18 * t13, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t17 * t10 + t21, 0, 0, 0, 0, 0, -t23 * t12, t23 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t9 * t21, g(1) * t2 - g(2) * t4 + t12 * t21;];
taug_reg = t5;
