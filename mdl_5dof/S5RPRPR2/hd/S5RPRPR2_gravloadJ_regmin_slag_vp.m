% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR2
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
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = qJ(1) + pkin(8);
t13 = qJ(3) + t15;
t10 = cos(t13);
t9 = sin(t13);
t22 = -t9 * pkin(3) + t10 * qJ(4);
t5 = g(2) * t9 - g(3) * t10;
t6 = g(2) * t10 + g(3) * t9;
t18 = sin(qJ(1));
t19 = cos(qJ(1));
t21 = g(2) * t19 + g(3) * t18;
t20 = -t10 * pkin(3) - t9 * qJ(4);
t14 = pkin(9) + qJ(5);
t12 = cos(t14);
t11 = sin(t14);
t4 = t6 * cos(pkin(9));
t3 = t6 * sin(pkin(9));
t2 = t6 * t12;
t1 = t6 * t11;
t7 = [0, t21, -g(2) * t18 + g(3) * t19, t21 * pkin(1), 0, t6, -t5, t4, -t3, t5, -g(2) * (-pkin(2) * cos(t15) - t19 * pkin(1) + t20) - g(3) * (-pkin(2) * sin(t15) - t18 * pkin(1) + t22), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t6, -t5, t4, -t3, t5, -g(2) * t20 - g(3) * t22, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 - t5 * t11, g(1) * t11 - t5 * t12;];
taug_reg = t7;
