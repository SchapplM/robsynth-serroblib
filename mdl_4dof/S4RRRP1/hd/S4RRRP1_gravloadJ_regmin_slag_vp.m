% Calculate minimal parameter regressor of gravitation load for
% S4RRRP1
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
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:24:20
% EndTime: 2019-05-04 19:24:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (70->15), mult. (42->23), div. (0->0), fcn. (36->6), ass. (0->15)
t12 = qJ(1) + qJ(2);
t10 = cos(t12);
t11 = qJ(3) + t12;
t8 = cos(t11);
t16 = pkin(2) * t10 + pkin(3) * t8;
t7 = sin(t11);
t9 = sin(t12);
t15 = -pkin(2) * t9 - pkin(3) * t7;
t1 = g(1) * t7 - g(2) * t8;
t14 = cos(qJ(1));
t13 = sin(qJ(1));
t4 = g(1) * t10 + g(2) * t9;
t3 = g(1) * t9 - g(2) * t10;
t2 = g(1) * t8 + g(2) * t7;
t5 = [0, g(1) * t13 - g(2) * t14, g(1) * t14 + g(2) * t13, 0, t3, t4, 0, t1, t2, -g(1) * (-pkin(1) * t13 + t15) - g(2) * (pkin(1) * t14 + t16); 0, 0, 0, 0, t3, t4, 0, t1, t2, -g(1) * t15 - g(2) * t16; 0, 0, 0, 0, 0, 0, 0, t1, t2, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
taug_reg  = t5;
