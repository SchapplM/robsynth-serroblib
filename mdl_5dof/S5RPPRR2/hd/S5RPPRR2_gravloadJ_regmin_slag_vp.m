% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = sin(qJ(1));
t18 = cos(qJ(1));
t19 = t18 * pkin(1) + t17 * qJ(2);
t14 = pkin(8) + qJ(4);
t4 = g(1) * t18 + g(2) * t17;
t3 = g(1) * t17 - g(2) * t18;
t11 = t18 * qJ(2);
t9 = qJ(5) + t14;
t8 = cos(t14);
t7 = sin(t14);
t6 = cos(t9);
t5 = sin(t9);
t2 = g(3) * t5 - t3 * t6;
t1 = g(3) * t6 + t3 * t5;
t10 = [0, t3, t4, -t3, -t4, -g(1) * (-t17 * pkin(1) + t11) - g(2) * t19, -t4 * sin(pkin(8)), -t4 * cos(pkin(8)), t3, -g(1) * (t11 + (-pkin(1) - qJ(3)) * t17) - g(2) * (t18 * qJ(3) + t19), 0, 0, 0, 0, 0, -t4 * t7, -t4 * t8, 0, 0, 0, 0, 0, -t4 * t5, -t4 * t6; 0, 0, 0, 0, 0, -t3, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t7 - t3 * t8, g(3) * t8 + t3 * t7, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t10;
