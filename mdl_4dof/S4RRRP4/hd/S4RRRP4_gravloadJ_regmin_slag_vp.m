% Calculate minimal parameter regressor of gravitation load for
% S4RRRP4
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
% taug_reg [4x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = cos(qJ(2));
t12 = qJ(2) + qJ(3);
t9 = cos(t12);
t17 = t15 * pkin(2) + pkin(3) * t9;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t6 = g(1) * t16 + g(2) * t14;
t5 = g(1) * t14 - g(2) * t16;
t8 = sin(t12);
t1 = -g(3) * t9 + t6 * t8;
t13 = sin(qJ(2));
t11 = -qJ(4) - pkin(6) - pkin(5);
t3 = pkin(1) + t17;
t2 = g(3) * t8 + t6 * t9;
t4 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t15, -t5 * t13, 0, 0, 0, 0, 0, t5 * t9, -t5 * t8, -t6, -g(1) * (-t16 * t11 - t14 * t3) - g(2) * (-t14 * t11 + t16 * t3); 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t15 + t6 * t13, g(3) * t13 + t6 * t15, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t17 - t6 * (-t13 * pkin(2) - pkin(3) * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t4;
