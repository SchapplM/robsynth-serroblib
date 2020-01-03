% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t6 = g(1) * t16 + g(2) * t13;
t25 = -g(3) * t12 + t6 * t15;
t24 = t16 * pkin(1) + t13 * qJ(2);
t22 = g(3) * t15;
t11 = sin(qJ(5));
t21 = t13 * t11;
t14 = cos(qJ(5));
t20 = t13 * t14;
t19 = t16 * t11;
t18 = t16 * t14;
t5 = g(1) * t13 - g(2) * t16;
t8 = t16 * qJ(2);
t4 = t12 * t18 - t21;
t3 = -t12 * t19 - t20;
t2 = -t12 * t20 - t19;
t1 = t12 * t21 - t18;
t7 = [0, t5, t6, -t5, -t6, -g(1) * (-t13 * pkin(1) + t8) - g(2) * t24, -t6, t5, -g(1) * (t8 + (-pkin(1) - qJ(3)) * t13) - g(2) * (t16 * qJ(3) + t24), 0, 0, 0, 0, 0, t5 * t12, t5 * t15, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, -t5, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t6 * t12 + t22, 0, 0, 0, 0, 0, -t25 * t14, t25 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t11 * t22, g(1) * t4 - g(2) * t2 + t14 * t22;];
taug_reg = t7;
