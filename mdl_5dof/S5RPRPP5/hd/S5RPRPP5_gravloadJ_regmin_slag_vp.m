% Calculate minimal parameter regressor of gravitation load for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = -pkin(1) - pkin(6);
t33 = -pkin(3) - pkin(4);
t21 = sin(qJ(1));
t20 = sin(qJ(3));
t26 = qJ(4) * t20;
t22 = cos(qJ(3));
t28 = t21 * t22;
t32 = pkin(3) * t28 + t21 * t26;
t23 = cos(qJ(1));
t31 = g(2) * t23;
t30 = t20 * t21;
t29 = t20 * t23;
t27 = t23 * pkin(1) + t21 * qJ(2);
t14 = t22 * qJ(4);
t25 = pkin(3) * t30 + t23 * pkin(6) + t27;
t6 = g(1) * t23 + g(2) * t21;
t15 = t23 * qJ(2);
t24 = pkin(3) * t29 - t23 * t14 + t15;
t5 = g(1) * t21 - t31;
t4 = t6 * t22;
t3 = t6 * t20;
t2 = g(1) * t28 - g(3) * t20 - t22 * t31;
t1 = g(1) * t30 - g(2) * t29 + g(3) * t22;
t7 = [0, t5, t6, -t5, -t6, -g(1) * (-t21 * pkin(1) + t15) - g(2) * t27, 0, 0, 0, 0, 0, -t3, -t4, -t3, t5, t4, -g(1) * (t34 * t21 + t24) - g(2) * (-t21 * t14 + t25), -t3, t4, -t5, -g(1) * (pkin(4) * t29 + t24) - g(2) * (-t23 * qJ(5) + t25) + (-g(1) * (qJ(5) + t34) - g(2) * (t20 * pkin(4) - t14)) * t21; 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, -g(1) * t32 - g(3) * (-t20 * pkin(3) + t14) - (-pkin(3) * t22 - t26) * t31, -t2, -t1, 0, -g(1) * (pkin(4) * t28 + t32) - g(3) * (t33 * t20 + t14) - (t33 * t22 - t26) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6;];
taug_reg = t7;
