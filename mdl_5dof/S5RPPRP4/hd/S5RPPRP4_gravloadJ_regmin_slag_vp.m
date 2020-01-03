% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = cos(qJ(1));
t27 = sin(qJ(1));
t26 = t28 * pkin(1) + t27 * qJ(2);
t25 = cos(pkin(7));
t24 = sin(pkin(7));
t23 = t28 * pkin(2) + t26;
t2 = -t27 * t24 - t28 * t25;
t3 = t28 * t24 - t27 * t25;
t22 = g(1) * t3 - g(2) * t2;
t21 = g(1) * t2 + g(2) * t3;
t20 = -t27 * pkin(1) + t28 * qJ(2);
t19 = -t27 * pkin(2) + t20;
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t18 = g(3) * t17 - t21 * t16;
t15 = -qJ(5) - pkin(6);
t9 = t17 * pkin(4) + pkin(3);
t5 = g(1) * t28 + g(2) * t27;
t4 = g(1) * t27 - g(2) * t28;
t1 = [0, t4, t5, t4, -t5, -g(1) * t20 - g(2) * t26, -t22, t21, -g(1) * t19 - g(2) * t23, 0, 0, 0, 0, 0, -t22 * t17, t22 * t16, -t21, -g(1) * (-t2 * t15 + t3 * t9 + t19) - g(2) * (-t3 * t15 - t2 * t9 + t23); 0, 0, 0, 0, 0, -t4, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -g(3) * t16 - t21 * t17, 0, t18 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22;];
taug_reg = t1;
