% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (68->33), mult. (108->42), div. (0->0), fcn. (112->8), ass. (0->22)
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t5 = g(1) * t17 - g(2) * t19;
t13 = pkin(8) + qJ(4);
t7 = sin(t13);
t8 = cos(t13);
t28 = -g(3) * t7 + t5 * t8;
t26 = g(3) * t8;
t25 = t19 * pkin(1) + t17 * qJ(2);
t16 = sin(qJ(5));
t24 = t17 * t16;
t18 = cos(qJ(5));
t23 = t17 * t18;
t22 = t19 * t16;
t21 = t19 * t18;
t6 = g(1) * t19 + g(2) * t17;
t10 = t19 * qJ(2);
t4 = t7 * t21 - t24;
t3 = t7 * t22 + t23;
t2 = t7 * t23 + t22;
t1 = -t7 * t24 + t21;
t9 = [0, t5, t6, -t5, -t6, -g(1) * (-t17 * pkin(1) + t10) - g(2) * t25, -t6 * sin(pkin(8)), -t6 * cos(pkin(8)), t5, -g(1) * (t10 + (-pkin(1) - qJ(3)) * t17) - g(2) * (t19 * qJ(3) + t25), 0, 0, 0, 0, 0, -t6 * t7, -t6 * t8, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t5 * t7 + t26, 0, 0, 0, 0, 0, -t28 * t18, t28 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t16 * t26, g(1) * t2 - g(2) * t4 + t18 * t26;];
taug_reg = t9;
