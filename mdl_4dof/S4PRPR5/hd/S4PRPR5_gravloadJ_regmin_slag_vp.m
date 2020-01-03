% Calculate minimal parameter regressor of gravitation load for
% S4PRPR5
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
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t3 = qJ(2) + pkin(7);
t1 = sin(t3);
t4 = sin(pkin(6));
t5 = cos(pkin(6));
t11 = g(1) * t5 + g(2) * t4;
t2 = cos(t3);
t18 = -g(3) * t2 + t11 * t1;
t17 = g(3) * t1;
t6 = sin(qJ(4));
t15 = t4 * t6;
t8 = cos(qJ(4));
t14 = t4 * t8;
t13 = t5 * t6;
t12 = t5 * t8;
t7 = sin(qJ(2));
t9 = cos(qJ(2));
t10 = -g(3) * t9 + t11 * t7;
t16 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t10, g(3) * t7 + t11 * t9, t10 * pkin(2), 0, 0, 0, 0, 0, t18 * t8, -t18 * t6; 0, 0, 0, 0, -g(1) * t4 + g(2) * t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t2 * t13 + t14) - g(2) * (-t2 * t15 - t12) + t6 * t17, -g(1) * (-t2 * t12 - t15) - g(2) * (-t2 * t14 + t13) + t8 * t17;];
taug_reg = t16;
