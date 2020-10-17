% Calculate minimal parameter regressor of gravitation load for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (47->21), mult. (85->32), div. (0->0), fcn. (86->8), ass. (0->14)
t7 = sin(pkin(6));
t9 = cos(pkin(6));
t14 = g(1) * t9 + g(2) * t7;
t10 = sin(qJ(2));
t11 = cos(qJ(2));
t1 = -g(3) * t11 + t14 * t10;
t18 = g(3) * t10;
t16 = t11 * t7;
t15 = t11 * t9;
t5 = pkin(7) + qJ(4);
t4 = cos(t5);
t3 = sin(t5);
t2 = t14 * t11 + t18;
t6 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t1, t2, t1 * cos(pkin(7)), -t1 * sin(pkin(7)), -t2, -g(3) * (t11 * pkin(2) + t10 * qJ(3)) + t14 * (pkin(2) * t10 - qJ(3) * t11), 0, 0, 0, 0, 0, t1 * t4, -t1 * t3; 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t3 * t15 + t7 * t4) - g(2) * (-t3 * t16 - t9 * t4) + t3 * t18, -g(1) * (-t4 * t15 - t7 * t3) - g(2) * (-t4 * t16 + t9 * t3) + t4 * t18;];
taug_reg = t6;
