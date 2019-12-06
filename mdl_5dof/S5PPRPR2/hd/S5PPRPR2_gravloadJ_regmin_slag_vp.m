% Calculate minimal parameter regressor of gravitation load for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t14 = g(1) * t9 + g(2) * t8;
t7 = pkin(8) + qJ(3);
t5 = sin(t7);
t6 = cos(t7);
t2 = g(3) * t5 + t14 * t6;
t19 = g(3) * t6;
t10 = sin(qJ(5));
t18 = t8 * t10;
t11 = cos(qJ(5));
t17 = t8 * t11;
t16 = t9 * t10;
t15 = t9 * t11;
t3 = -g(1) * t8 + g(2) * t9;
t1 = t14 * t5 - t19;
t4 = [-g(3), -g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, t3, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t1, t2, -t1, -t2, -g(3) * (t6 * pkin(3) + t5 * qJ(4)) + t14 * (pkin(3) * t5 - qJ(4) * t6), 0, 0, 0, 0, 0, -t2 * t10, -t2 * t11; 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t5 * t15 - t18) - g(2) * (t5 * t17 + t16) + t11 * t19, -g(1) * (-t5 * t16 - t17) - g(2) * (-t5 * t18 + t15) - t10 * t19;];
taug_reg = t4;
