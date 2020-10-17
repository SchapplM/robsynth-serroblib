% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (95->23), mult. (96->31), div. (0->0), fcn. (112->8), ass. (0->18)
t16 = qJ(1) + pkin(8);
t14 = sin(t16);
t15 = cos(t16);
t17 = sin(qJ(4));
t18 = cos(qJ(4));
t1 = -t14 * t17 - t15 * t18;
t2 = -t14 * t18 + t15 * t17;
t13 = g(1) * t2 - g(2) * t1;
t7 = sin(qJ(5));
t20 = t13 * t7;
t9 = cos(qJ(5));
t19 = t13 * t9;
t12 = g(1) * t1 + g(2) * t2;
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t11 = g(1) * t8 - g(2) * t10;
t3 = g(1) * t14 - g(2) * t15;
t4 = [0, t11, g(1) * t10 + g(2) * t8, t11 * pkin(1), t3, -g(1) * t15 - g(2) * t14, -g(1) * (-t8 * pkin(1) - t14 * pkin(2) + t15 * qJ(3)) - g(2) * (t10 * pkin(1) + t15 * pkin(2) + t14 * qJ(3)), 0, -t13, t12, 0, 0, 0, 0, 0, -t19, t20; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, 0, 0, 0, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t9 - t12 * t7, -g(3) * t7 - t12 * t9;];
taug_reg = t4;
