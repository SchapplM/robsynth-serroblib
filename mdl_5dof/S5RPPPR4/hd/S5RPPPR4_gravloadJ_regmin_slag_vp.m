% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (76->23), mult. (64->29), div. (0->0), fcn. (58->8), ass. (0->14)
t13 = qJ(1) + pkin(7);
t10 = cos(t13);
t17 = cos(qJ(1));
t8 = sin(t13);
t20 = t17 * pkin(1) + t10 * pkin(2) + t8 * qJ(3);
t16 = sin(qJ(1));
t19 = -t16 * pkin(1) + t10 * qJ(3);
t1 = g(1) * t8 - g(2) * t10;
t2 = -g(1) * t10 - g(2) * t8;
t18 = g(1) * t16 - g(2) * t17;
t12 = pkin(8) + qJ(5);
t9 = cos(t12);
t7 = sin(t12);
t3 = [0, t18, g(1) * t17 + g(2) * t16, t18 * pkin(1), -t1, t2, -g(1) * (-t8 * pkin(2) + t19) - g(2) * t20, t2 * sin(pkin(8)), t2 * cos(pkin(8)), t1, -g(1) * ((-pkin(2) - qJ(4)) * t8 + t19) - g(2) * (t10 * qJ(4) + t20), 0, 0, 0, 0, 0, t2 * t7, t2 * t9; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t7 - t1 * t9, g(3) * t9 + t1 * t7;];
taug_reg = t3;
