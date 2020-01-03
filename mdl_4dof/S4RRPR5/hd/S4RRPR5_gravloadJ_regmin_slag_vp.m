% Calculate minimal parameter regressor of gravitation load for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% taug_reg [4x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t11 = qJ(1) + qJ(2);
t10 = cos(t11);
t9 = sin(t11);
t17 = t10 * pkin(2) + t9 * qJ(3);
t16 = -t9 * pkin(2) + t10 * qJ(3);
t3 = g(1) * t9 - g(2) * t10;
t4 = g(1) * t10 + g(2) * t9;
t15 = cos(qJ(1));
t14 = cos(qJ(4));
t13 = sin(qJ(1));
t12 = sin(qJ(4));
t2 = t4 * t14;
t1 = t4 * t12;
t5 = [0, g(1) * t13 - g(2) * t15, g(1) * t15 + g(2) * t13, 0, t3, t4, -t3, -t4, -g(1) * (-t13 * pkin(1) + t16) - g(2) * (t15 * pkin(1) + t17), 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, t3, t4, -t3, -t4, -g(1) * t16 - g(2) * t17, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t12 - t3 * t14, g(3) * t14 + t3 * t12;];
taug_reg = t5;
