% Calculate minimal parameter regressor of gravitation load for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(1));
t17 = cos(qJ(1));
t23 = t17 * pkin(1) + t16 * qJ(2);
t22 = pkin(6) + qJ(4);
t21 = cos(t22);
t20 = sin(t22);
t1 = -t16 * t20 - t17 * t21;
t2 = -t16 * t21 + t17 * t20;
t19 = g(1) * t2 - g(2) * t1;
t18 = g(1) * t1 + g(2) * t2;
t15 = cos(pkin(6));
t14 = sin(pkin(6));
t11 = t17 * qJ(2);
t6 = g(1) * t17 + g(2) * t16;
t5 = g(1) * t16 - g(2) * t17;
t4 = t16 * t14 + t17 * t15;
t3 = t17 * t14 - t16 * t15;
t7 = [0, t5, t6, t5, -t6, -g(1) * (-t16 * pkin(1) + t11) - g(2) * t23, -g(1) * t3 - g(2) * t4, -g(1) * t4 + g(2) * t3, -g(1) * (t11 + (-pkin(1) - pkin(2)) * t16) - g(2) * (t17 * pkin(2) + t23) 0, -t19, t18; 0, 0, 0, 0, 0, -t5, 0, 0, -t5, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18;];
taug_reg  = t7;
