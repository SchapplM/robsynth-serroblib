% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = cos(qJ(1));
t27 = sin(qJ(1));
t26 = t28 * pkin(1) + t27 * qJ(2);
t25 = cos(pkin(7));
t24 = sin(pkin(7));
t23 = t28 * pkin(2) + t26;
t2 = -t27 * t24 - t28 * t25;
t3 = t28 * t24 - t27 * t25;
t22 = g(1) * t3 - g(2) * t2;
t21 = g(1) * t2 + g(2) * t3;
t20 = -t27 * pkin(1) + t28 * qJ(2);
t19 = -t27 * pkin(2) + t20;
t16 = pkin(8) + qJ(5);
t10 = cos(t16);
t9 = sin(t16);
t5 = g(1) * t28 + g(2) * t27;
t4 = g(1) * t27 - g(2) * t28;
t1 = [0, t4, t5, t4, -t5, -g(1) * t20 - g(2) * t26, -t22, t21, -g(1) * t19 - g(2) * t23, -t22 * cos(pkin(8)), t22 * sin(pkin(8)), -t21, -g(1) * (t3 * pkin(3) + t2 * qJ(4) + t19) - g(2) * (-t2 * pkin(3) + t3 * qJ(4) + t23), 0, 0, 0, 0, 0, -t22 * t10, t22 * t9; 0, 0, 0, 0, 0, -t4, 0, 0, -t4, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t10 - t21 * t9, -g(3) * t9 - t21 * t10;];
taug_reg = t1;
