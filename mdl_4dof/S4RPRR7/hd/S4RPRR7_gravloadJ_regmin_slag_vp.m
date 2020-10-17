% Calculate minimal parameter regressor of gravitation load for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (59->23), mult. (92->38), div. (0->0), fcn. (98->8), ass. (0->20)
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t6 = g(1) * t15 + g(2) * t13;
t9 = pkin(7) + qJ(3);
t7 = sin(t9);
t8 = cos(t9);
t16 = -g(3) * t8 + t6 * t7;
t22 = g(3) * t7;
t12 = sin(qJ(4));
t20 = t13 * t12;
t14 = cos(qJ(4));
t19 = t13 * t14;
t18 = t15 * t12;
t17 = t15 * t14;
t5 = g(1) * t13 - g(2) * t15;
t4 = t8 * t17 + t20;
t3 = -t8 * t18 + t19;
t2 = -t8 * t19 + t18;
t1 = t8 * t20 + t17;
t10 = [0, t5, t6, t5 * cos(pkin(7)), -t5 * sin(pkin(7)), -t6, -g(1) * (-t13 * pkin(1) + t15 * qJ(2)) - g(2) * (t15 * pkin(1) + t13 * qJ(2)), 0, 0, 0, 0, 0, t5 * t8, -t5 * t7, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t6 * t8 + t22, 0, 0, 0, 0, 0, t16 * t14, -t16 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t12 * t22, g(1) * t4 - g(2) * t2 + t14 * t22;];
taug_reg = t10;
