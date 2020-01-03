% Calculate minimal parameter regressor of gravitation load for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% taug_reg [4x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t3 = g(1) * t15 + g(2) * t13;
t10 = qJ(2) + pkin(6);
t6 = sin(t10);
t7 = cos(t10);
t18 = t7 * pkin(3) + t6 * qJ(4);
t2 = g(1) * t13 - g(2) * t15;
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t16 = -g(3) * t14 + t3 * t12;
t11 = -qJ(3) - pkin(5);
t8 = t14 * pkin(2);
t5 = t8 + pkin(1);
t4 = t15 * t5;
t1 = -g(3) * t7 + t3 * t6;
t9 = [0, t2, t3, 0, 0, 0, 0, 0, t2 * t14, -t2 * t12, -t3, -g(1) * (-t15 * t11 - t13 * t5) - g(2) * (-t13 * t11 + t4), t2 * t7, -t3, t2 * t6, -g(2) * t4 + (g(1) * t11 - g(2) * t18) * t15 + (-g(1) * (-t18 - t5) + g(2) * t11) * t13; 0, 0, 0, 0, 0, 0, 0, 0, t16, g(3) * t12 + t3 * t14, 0, t16 * pkin(2), t1, 0, -g(3) * t6 - t3 * t7, -g(3) * (t18 + t8) + t3 * (pkin(2) * t12 + pkin(3) * t6 - qJ(4) * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
