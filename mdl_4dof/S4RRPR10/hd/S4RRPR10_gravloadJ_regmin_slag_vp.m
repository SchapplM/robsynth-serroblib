% Calculate minimal parameter regressor of gravitation load for
% S4RRPR10
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t9 = g(1) * t16 + g(2) * t13;
t12 = sin(qJ(2));
t15 = cos(qJ(2));
t2 = g(3) * t12 + t9 * t15;
t26 = g(3) * t15;
t11 = sin(qJ(4));
t25 = t13 * t11;
t14 = cos(qJ(4));
t24 = t13 * t14;
t23 = t16 * t11;
t22 = t16 * t14;
t21 = g(1) * t13 - g(2) * t16;
t20 = t15 * pkin(2) + t12 * qJ(3);
t18 = pkin(1) + t20;
t8 = t21 * t15;
t7 = t21 * t12;
t6 = -t12 * t25 + t22;
t5 = t12 * t24 + t23;
t4 = t12 * t23 + t24;
t3 = t12 * t22 - t25;
t1 = t9 * t12 - t26;
t10 = [0, t21, t9, 0, 0, 0, 0, 0, t8, -t7, -t9, -t8, t7, (-g(1) * pkin(5) - g(2) * t18) * t16 + (-g(2) * pkin(5) + g(1) * t18) * t13, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(3) * t20 + t9 * (pkin(2) * t12 - qJ(3) * t15), 0, 0, 0, 0, 0, -t2 * t11, -t2 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t14 * t26, g(1) * t4 - g(2) * t6 - t11 * t26;];
taug_reg = t10;
