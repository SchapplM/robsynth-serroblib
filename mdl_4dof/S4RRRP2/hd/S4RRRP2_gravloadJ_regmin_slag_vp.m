% Calculate minimal parameter regressor of gravitation load for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:09
% DurationCPUTime: 0.09s
% Computational Cost: add. (71->20), mult. (70->26), div. (0->0), fcn. (63->6), ass. (0->17)
t10 = -qJ(4) - pkin(6);
t13 = cos(qJ(3));
t6 = t13 * pkin(3) + pkin(2);
t9 = qJ(1) + qJ(2);
t7 = sin(t9);
t8 = cos(t9);
t17 = -t7 * t10 + t8 * t6;
t4 = g(1) * t8 + g(2) * t7;
t3 = g(1) * t7 - g(2) * t8;
t16 = -t8 * t10 - t7 * t6;
t11 = sin(qJ(3));
t15 = -g(3) * t13 + t4 * t11;
t14 = cos(qJ(1));
t12 = sin(qJ(1));
t2 = t3 * t13;
t1 = t3 * t11;
t5 = [0, g(1) * t12 - g(2) * t14, g(1) * t14 + g(2) * t12, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t12 * pkin(1) + t16) - g(2) * (t14 * pkin(1) + t17); 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t16 - g(2) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, g(3) * t11 + t4 * t13, 0, t15 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg = t5;
