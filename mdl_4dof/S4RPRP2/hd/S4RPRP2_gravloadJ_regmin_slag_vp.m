% Calculate minimal parameter regressor of gravitation load for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:15:22
% EndTime: 2019-05-04 19:15:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (33->21), mult. (66->26), div. (0->0), fcn. (68->4), ass. (0->16)
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t21 = t16 * pkin(1) + t14 * qJ(2);
t13 = sin(qJ(3));
t20 = t13 * t16;
t19 = t14 * t13;
t15 = cos(qJ(3));
t1 = -t16 * t15 - t19;
t2 = -t14 * t15 + t20;
t18 = g(1) * t2 - g(2) * t1;
t17 = g(1) * t1 + g(2) * t2;
t10 = t16 * qJ(2);
t8 = pkin(3) * t15 + pkin(2);
t4 = g(1) * t16 + g(2) * t14;
t3 = g(1) * t14 - g(2) * t16;
t5 = [0, t3, t4, t3, -t4, -g(1) * (-t14 * pkin(1) + t10) - g(2) * t21, 0, -t18, t17, -g(1) * (pkin(3) * t20 + t10 + (-pkin(1) - t8) * t14) - g(2) * (pkin(3) * t19 + t16 * t8 + t21); 0, 0, 0, 0, 0, -t3, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, t18, -t17, t18 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3);];
taug_reg  = t5;
