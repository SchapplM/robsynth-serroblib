% Calculate minimal parameter regressor of gravitation load for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (38->24), mult. (99->38), div. (0->0), fcn. (105->6), ass. (0->21)
t16 = sin(qJ(1));
t23 = g(1) * t16;
t18 = cos(qJ(1));
t22 = t18 * pkin(1) + t16 * qJ(2);
t8 = g(1) * t18 + g(2) * t16;
t7 = -g(2) * t18 + t23;
t13 = sin(pkin(6));
t14 = cos(pkin(6));
t21 = pkin(2) * t14 + qJ(3) * t13;
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t20 = t13 * t17 - t14 * t15;
t19 = t13 * t15 + t14 * t17;
t10 = t18 * qJ(2);
t6 = t7 * t14;
t5 = t7 * t13;
t4 = t19 * t18;
t3 = t20 * t18;
t2 = t19 * t16;
t1 = t20 * t16;
t9 = [0, t7, t8, t6, -t5, -t8, -g(1) * (-t16 * pkin(1) + t10) - g(2) * t22, t6, -t8, t5, -g(1) * t10 - g(2) * (t21 * t18 + t22) - (-pkin(1) - t21) * t23, 0, 0, 0, 0, 0, g(1) * t2 - g(2) * t4, g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t14 - t8 * t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t1 + g(3) * t19, g(1) * t4 + g(2) * t2 + g(3) * t20;];
taug_reg = t9;
