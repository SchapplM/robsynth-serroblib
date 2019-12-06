% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t19 = t18 * pkin(1) + t16 * qJ(2);
t6 = g(1) * t18 + g(2) * t16;
t5 = g(1) * t16 - g(2) * t18;
t14 = qJ(3) + qJ(4);
t7 = sin(t14);
t8 = cos(t14);
t2 = g(3) * t7 - t5 * t8;
t17 = cos(qJ(3));
t15 = sin(qJ(3));
t13 = -qJ(5) - pkin(7) - pkin(6);
t10 = t18 * qJ(2);
t3 = t15 * pkin(3) + pkin(4) * t7;
t1 = g(3) * t8 + t5 * t7;
t4 = [0, t5, t6, -t5, -t6, -g(1) * (-t16 * pkin(1) + t10) - g(2) * t19, 0, 0, 0, 0, 0, -t6 * t15, -t6 * t17, 0, 0, 0, 0, 0, -t6 * t7, -t6 * t8, t5, -g(1) * (t18 * t3 + t10 + (-pkin(1) + t13) * t16) - g(2) * (-t18 * t13 + t16 * t3 + t19); 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t15 - t5 * t17, g(3) * t17 + t5 * t15, 0, 0, 0, 0, 0, t2, t1, 0, g(3) * t3 - t5 * (t17 * pkin(3) + pkin(4) * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t4;
