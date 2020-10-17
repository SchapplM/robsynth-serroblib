% Calculate minimal parameter regressor of gravitation load for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% taug_reg [4x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (32->24), mult. (86->36), div. (0->0), fcn. (92->6), ass. (0->19)
t12 = cos(qJ(3));
t10 = sin(qJ(1));
t13 = cos(qJ(1));
t5 = g(1) * t10 - g(2) * t13;
t9 = sin(qJ(3));
t21 = -g(3) * t9 + t5 * t12;
t19 = g(3) * t12;
t8 = sin(qJ(4));
t18 = t10 * t8;
t17 = t13 * t8;
t11 = cos(qJ(4));
t16 = t10 * t11;
t15 = t13 * t11;
t6 = g(1) * t13 + g(2) * t10;
t4 = t9 * t15 - t18;
t3 = t9 * t17 + t16;
t2 = t9 * t16 + t17;
t1 = -t9 * t18 + t15;
t7 = [0, t5, t6, -t5, -t6, -g(1) * (-t10 * pkin(1) + t13 * qJ(2)) - g(2) * (t13 * pkin(1) + t10 * qJ(2)), 0, 0, 0, 0, 0, -t6 * t9, -t6 * t12, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t5 * t9 + t19, 0, 0, 0, 0, 0, -t21 * t11, t21 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t8 * t19, g(1) * t2 - g(2) * t4 + t11 * t19;];
taug_reg = t7;
