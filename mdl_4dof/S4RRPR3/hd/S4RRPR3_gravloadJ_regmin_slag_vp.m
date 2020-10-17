% Calculate minimal parameter regressor of gravitation load for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (55->16), mult. (48->25), div. (0->0), fcn. (44->8), ass. (0->17)
t10 = qJ(1) + qJ(2);
t7 = pkin(7) + t10;
t5 = sin(t7);
t6 = cos(t7);
t16 = g(1) * t6 + g(2) * t5;
t15 = g(1) * t5 - g(2) * t6;
t8 = sin(t10);
t9 = cos(t10);
t3 = g(1) * t8 - g(2) * t9;
t14 = cos(qJ(1));
t13 = cos(qJ(4));
t12 = sin(qJ(1));
t11 = sin(qJ(4));
t4 = g(1) * t9 + g(2) * t8;
t2 = t15 * t13;
t1 = t15 * t11;
t17 = [0, g(1) * t12 - g(2) * t14, g(1) * t14 + g(2) * t12, 0, t3, t4, -g(1) * (-t12 * pkin(1) - pkin(2) * t8) - g(2) * (t14 * pkin(1) + pkin(2) * t9), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t3, t4, t3 * pkin(2), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t13 + t16 * t11, g(3) * t11 + t16 * t13;];
taug_reg = t17;
