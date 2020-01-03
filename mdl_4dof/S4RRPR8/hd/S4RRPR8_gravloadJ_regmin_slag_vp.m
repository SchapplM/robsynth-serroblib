% Calculate minimal parameter regressor of gravitation load for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = sin(qJ(4));
t14 = sin(qJ(2));
t16 = cos(qJ(4));
t17 = cos(qJ(2));
t22 = t17 * t13 - t14 * t16;
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t11 = g(1) * t18 + g(2) * t15;
t25 = g(1) * t15 - g(2) * t18;
t24 = t17 * pkin(2) + t14 * qJ(3);
t9 = t14 * t13 + t17 * t16;
t21 = pkin(1) + t24;
t3 = t22 * t15;
t5 = t22 * t18;
t20 = g(1) * t5 + g(2) * t3 + g(3) * t9;
t4 = t9 * t15;
t6 = t9 * t18;
t19 = g(1) * t6 + g(2) * t4 - g(3) * t22;
t8 = t25 * t17;
t7 = t25 * t14;
t2 = g(3) * t14 + t11 * t17;
t1 = -g(3) * t17 + t11 * t14;
t10 = [0, t25, t11, 0, 0, 0, 0, 0, t8, -t7, t8, -t11, t7, (-g(1) * pkin(5) - g(2) * t21) * t18 + (-g(2) * pkin(5) + g(1) * t21) * t15, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 + g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t24 + t11 * (pkin(2) * t14 - qJ(3) * t17), 0, 0, 0, 0, 0, -t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t19;];
taug_reg = t10;
