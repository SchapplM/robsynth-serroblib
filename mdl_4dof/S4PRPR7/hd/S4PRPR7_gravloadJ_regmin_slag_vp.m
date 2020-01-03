% Calculate minimal parameter regressor of gravitation load for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% taug_reg [4x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t4 = sin(pkin(6));
t5 = cos(pkin(6));
t12 = g(1) * t5 + g(2) * t4;
t7 = sin(qJ(2));
t9 = cos(qJ(2));
t2 = g(3) * t7 + t12 * t9;
t15 = g(3) * t9;
t6 = sin(qJ(4));
t14 = t6 * t7;
t8 = cos(qJ(4));
t13 = t7 * t8;
t1 = t12 * t7 - t15;
t3 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t1, t2, -t1, -t2, -g(3) * (t9 * pkin(2) + t7 * qJ(3)) + t12 * (pkin(2) * t7 - qJ(3) * t9), 0, 0, 0, 0, 0, -t2 * t6, -t2 * t8; 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t5 * t13 - t4 * t6) - g(2) * (t4 * t13 + t5 * t6) + t8 * t15, -g(1) * (-t5 * t14 - t4 * t8) - g(2) * (-t4 * t14 + t5 * t8) - t6 * t15;];
taug_reg = t3;
