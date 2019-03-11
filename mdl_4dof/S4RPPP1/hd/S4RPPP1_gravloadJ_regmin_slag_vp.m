% Calculate minimal parameter regressor of gravitation load for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(pkin(4));
t34 = pkin(3) * t17;
t33 = g(3) * t17;
t16 = sin(pkin(6));
t20 = sin(qJ(1));
t32 = t20 * t16;
t18 = cos(pkin(6));
t31 = t20 * t18;
t21 = cos(qJ(1));
t30 = t21 * t16;
t29 = t21 * t18;
t27 = qJ(2) * t17;
t28 = t21 * pkin(1) + t20 * t27;
t26 = -t20 * pkin(1) + t21 * t27;
t19 = cos(pkin(4));
t7 = -t19 * t29 + t32;
t9 = t19 * t31 + t30;
t2 = g(1) * t7 - g(2) * t9;
t10 = -t19 * t32 + t29;
t8 = t19 * t30 + t31;
t3 = g(1) * t8 - g(2) * t10;
t25 = g(1) * t21 + g(2) * t20;
t24 = g(1) * t20 - g(2) * t21;
t23 = t10 * pkin(2) + t9 * qJ(3) + t28;
t22 = -t8 * pkin(2) - t7 * qJ(3) + t26;
t11 = t25 * t17;
t4 = -g(3) * t19 - t24 * t17;
t1 = -g(1) * t9 - g(2) * t7 + t18 * t33;
t5 = [0, t24, t25, t3, -t2, -t11, -g(1) * t26 - g(2) * t28, -t11, -t3, t2, -g(1) * t22 - g(2) * t23, -t11, t2, t3, -g(1) * (-t8 * qJ(4) + t21 * t34 + t22) - g(2) * (t10 * qJ(4) + t20 * t34 + t23); 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, t4, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8 - t16 * t33;];
taug_reg  = t5;
