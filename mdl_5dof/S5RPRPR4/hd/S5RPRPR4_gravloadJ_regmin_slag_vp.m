% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t9 = qJ(1) + pkin(8);
t6 = sin(t9);
t7 = cos(t9);
t18 = g(2) * t7 + g(3) * t6;
t17 = g(2) * t6 - g(3) * t7;
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t16 = -g(2) * t14 - g(3) * t12;
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t15 = -g(1) * t13 + t17 * t11;
t10 = -qJ(4) - pkin(6);
t8 = qJ(3) + pkin(9) + qJ(5);
t5 = t13 * pkin(3) + pkin(2);
t4 = cos(t8);
t3 = sin(t8);
t2 = -g(1) * t4 + t17 * t3;
t1 = g(1) * t3 + t17 * t4;
t19 = [0, t16, g(2) * t12 - g(3) * t14, t16 * pkin(1), 0, 0, 0, 0, 0, -t18 * t13, t18 * t11, -t17, -g(2) * (t14 * pkin(1) - t6 * t10 + t7 * t5) - g(3) * (t12 * pkin(1) + t7 * t10 + t6 * t5), 0, 0, 0, 0, 0, -t18 * t4, t18 * t3; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, g(1) * t11 + t17 * t13, 0, t15 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t19;
