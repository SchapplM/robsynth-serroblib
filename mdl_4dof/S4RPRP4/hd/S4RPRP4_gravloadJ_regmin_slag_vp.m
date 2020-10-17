% Calculate minimal parameter regressor of gravitation load for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (66->21), mult. (80->27), div. (0->0), fcn. (73->6), ass. (0->18)
t7 = qJ(1) + pkin(6);
t5 = sin(t7);
t6 = cos(t7);
t18 = g(1) * t6 + g(2) * t5;
t17 = g(1) * t5 - g(2) * t6;
t11 = cos(qJ(1));
t9 = sin(qJ(1));
t16 = g(1) * t9 - g(2) * t11;
t10 = cos(qJ(3));
t8 = sin(qJ(3));
t14 = t10 * pkin(3) + t8 * qJ(4);
t13 = pkin(2) + t14;
t12 = t16 * pkin(1);
t4 = t17 * t10;
t3 = t17 * t8;
t2 = g(3) * t8 + t18 * t10;
t1 = -g(3) * t10 + t18 * t8;
t15 = [0, t16, g(1) * t11 + g(2) * t9, t12, 0, 0, 0, 0, 0, t4, -t3, t4, -t18, t3, t12 + (-g(1) * pkin(5) - g(2) * t13) * t6 + (-g(2) * pkin(5) + g(1) * t13) * t5; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t14 + t18 * (pkin(3) * t8 - qJ(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t15;
