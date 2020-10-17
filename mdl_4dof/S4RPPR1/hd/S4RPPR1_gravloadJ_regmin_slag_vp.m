% Calculate minimal parameter regressor of gravitation load for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:08:42
% EndTime: 2019-05-04 19:08:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->19), mult. (46->25), div. (0->0), fcn. (50->6), ass. (0->14)
t16 = cos(qJ(4));
t15 = sin(qJ(4));
t9 = qJ(1) + pkin(6);
t7 = sin(t9);
t8 = cos(t9);
t1 = -t7 * t15 - t8 * t16;
t2 = t8 * t15 - t7 * t16;
t14 = g(1) * t2 - g(2) * t1;
t13 = g(1) * t1 + g(2) * t2;
t10 = sin(qJ(1));
t11 = cos(qJ(1));
t12 = g(1) * t10 - g(2) * t11;
t3 = g(1) * t7 - g(2) * t8;
t4 = [0, t12, g(1) * t11 + g(2) * t10, t12 * pkin(1), t3, -g(1) * t8 - g(2) * t7, -g(1) * (-pkin(1) * t10 - t7 * pkin(2) + t8 * qJ(3)) - g(2) * (pkin(1) * t11 + t8 * pkin(2) + t7 * qJ(3)) 0, -t14, t13; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0; 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13;];
taug_reg  = t4;
