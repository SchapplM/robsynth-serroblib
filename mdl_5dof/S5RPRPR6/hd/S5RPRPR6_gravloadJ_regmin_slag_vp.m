% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:49
% DurationCPUTime: 0.09s
% Computational Cost: add. (110->27), mult. (68->27), div. (0->0), fcn. (62->8), ass. (0->16)
t12 = qJ(1) + pkin(8);
t11 = qJ(3) + t12;
t10 = cos(t11);
t9 = sin(t11);
t19 = t10 * pkin(3) + t9 * qJ(4);
t18 = -t9 * pkin(3) + t10 * qJ(4);
t3 = g(1) * t9 - g(2) * t10;
t4 = g(1) * t10 + g(2) * t9;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t17 = g(1) * t14 - g(2) * t16;
t15 = cos(qJ(5));
t13 = sin(qJ(5));
t2 = t4 * t15;
t1 = t4 * t13;
t5 = [0, t17, g(1) * t16 + g(2) * t14, t17 * pkin(1), 0, t3, t4, -t3, -t4, -g(1) * (-pkin(2) * sin(t12) - t14 * pkin(1) + t18) - g(2) * (pkin(2) * cos(t12) + t16 * pkin(1) + t19), 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t3, t4, -t3, -t4, -g(1) * t18 - g(2) * t19, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t13 - t3 * t15, g(3) * t15 + t3 * t13;];
taug_reg = t5;
