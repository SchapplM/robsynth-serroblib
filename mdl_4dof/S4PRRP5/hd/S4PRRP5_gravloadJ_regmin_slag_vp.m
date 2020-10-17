% Calculate minimal parameter regressor of gravitation load for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x13]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (37->21), mult. (87->32), div. (0->0), fcn. (85->6), ass. (0->16)
t4 = sin(pkin(6));
t5 = cos(pkin(6));
t13 = g(1) * t5 + g(2) * t4;
t10 = cos(qJ(2));
t8 = sin(qJ(2));
t1 = -g(3) * t10 + t13 * t8;
t17 = g(3) * t8;
t7 = sin(qJ(3));
t15 = t10 * t7;
t9 = cos(qJ(3));
t14 = t10 * t9;
t11 = -g(1) * (-t5 * t15 + t4 * t9) - g(2) * (-t4 * t15 - t5 * t9) + t7 * t17;
t6 = -qJ(4) - pkin(5);
t3 = t9 * pkin(3) + pkin(2);
t2 = t13 * t10 + t17;
t12 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t9, -t1 * t7, -t2, -g(3) * (t10 * t3 - t8 * t6) + t13 * (t10 * t6 + t3 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -g(1) * (-t5 * t14 - t4 * t7) - g(2) * (-t4 * t14 + t5 * t7) + t9 * t17, 0, t11 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
