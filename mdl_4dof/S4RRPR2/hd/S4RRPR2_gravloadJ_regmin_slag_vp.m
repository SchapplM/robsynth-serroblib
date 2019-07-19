% Calculate minimal parameter regressor of gravitation load for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = cos(qJ(4));
t21 = sin(qJ(4));
t15 = qJ(1) + qJ(2);
t13 = sin(t15);
t14 = cos(t15);
t20 = t14 * pkin(2) + t13 * qJ(3);
t19 = -t13 * pkin(2) + t14 * qJ(3);
t3 = -t13 * t21 - t14 * t22;
t4 = -t13 * t22 + t14 * t21;
t18 = g(1) * t4 - g(2) * t3;
t2 = g(1) * t3 + g(2) * t4;
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t6 = g(1) * t14 + g(2) * t13;
t5 = g(1) * t13 - g(2) * t14;
t1 = [0, g(1) * t16 - g(2) * t17, g(1) * t17 + g(2) * t16, 0, t5, t6, t5, -t6, -g(1) * (-t16 * pkin(1) + t19) - g(2) * (t17 * pkin(1) + t20), 0, -t18, t2; 0, 0, 0, 0, t5, t6, t5, -t6, -g(1) * t19 - g(2) * t20, 0, -t18, t2; 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t2;];
taug_reg  = t1;
