% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [5x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (84->31), mult. (130->42), div. (0->0), fcn. (136->8), ass. (0->25)
t14 = qJ(3) + qJ(4);
t11 = sin(t14);
t12 = cos(t14);
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t9 = g(1) * t17 - g(2) * t20;
t27 = -g(3) * t11 + t9 * t12;
t25 = g(3) * t12;
t15 = sin(qJ(5));
t24 = t17 * t15;
t18 = cos(qJ(5));
t23 = t17 * t18;
t22 = t20 * t15;
t21 = t20 * t18;
t10 = g(1) * t20 + g(2) * t17;
t19 = cos(qJ(3));
t16 = sin(qJ(3));
t8 = t11 * t21 - t24;
t7 = t11 * t22 + t23;
t6 = t11 * t23 + t22;
t5 = -t11 * t24 + t21;
t3 = t9 * t11 + t25;
t2 = t27 * t18;
t1 = t27 * t15;
t4 = [0, t9, t10, -t9, -t10, -g(1) * (-t17 * pkin(1) + t20 * qJ(2)) - g(2) * (t20 * pkin(1) + t17 * qJ(2)), 0, 0, 0, 0, 0, -t10 * t16, -t10 * t19, 0, 0, 0, 0, 0, -t10 * t11, -t10 * t12, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t16 - t9 * t19, g(3) * t19 + t9 * t16, 0, 0, 0, 0, 0, -t27, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t15 * t25, g(1) * t6 - g(2) * t8 + t18 * t25;];
taug_reg = t4;
