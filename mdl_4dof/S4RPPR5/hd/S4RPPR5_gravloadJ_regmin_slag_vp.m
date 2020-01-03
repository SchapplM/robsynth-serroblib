% Calculate minimal parameter regressor of gravitation load for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% taug_reg [4x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t18 = t14 * pkin(1) + t12 * qJ(2);
t17 = cos(pkin(6));
t10 = sin(pkin(6));
t1 = t14 * t10 - t12 * t17;
t2 = t12 * t10 + t14 * t17;
t16 = g(1) * t2 - g(2) * t1;
t15 = g(1) * t1 + g(2) * t2;
t13 = cos(qJ(4));
t11 = sin(qJ(4));
t7 = t14 * qJ(2);
t4 = g(1) * t14 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t14;
t5 = [0, t3, t4, t3, -t4, -g(1) * (-t12 * pkin(1) + t7) - g(2) * t18, -t15, -t16, -g(1) * (t7 + (-pkin(1) - pkin(2)) * t12) - g(2) * (t14 * pkin(2) + t18), 0, 0, 0, 0, 0, -t15 * t13, t15 * t11; 0, 0, 0, 0, 0, -t3, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t13 + t16 * t11, -g(3) * t11 + t16 * t13;];
taug_reg = t5;
