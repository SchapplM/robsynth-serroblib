% Calculate minimal parameter regressor of gravitation load for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = cos(qJ(3));
t7 = qJ(1) + pkin(7);
t5 = sin(t7);
t6 = cos(t7);
t17 = g(1) * t6 + g(2) * t5;
t9 = sin(qJ(3));
t14 = -g(3) * t12 + t17 * t9;
t21 = g(3) * t9;
t8 = sin(qJ(4));
t19 = t12 * t8;
t11 = cos(qJ(4));
t18 = t11 * t12;
t16 = g(1) * t5 - g(2) * t6;
t10 = sin(qJ(1));
t13 = cos(qJ(1));
t15 = g(1) * t10 - g(2) * t13;
t4 = t6 * t18 + t5 * t8;
t3 = t5 * t11 - t6 * t19;
t2 = -t5 * t18 + t6 * t8;
t1 = t6 * t11 + t5 * t19;
t20 = [0, t15, g(1) * t13 + g(2) * t10, t15 * pkin(1), 0, 0, 0, 0, 0, t16 * t12, -t16 * t9, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t17 * t12 + t21, 0, 0, 0, 0, 0, t14 * t11, -t14 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t8 * t21, g(1) * t4 - g(2) * t2 + t11 * t21;];
taug_reg = t20;
