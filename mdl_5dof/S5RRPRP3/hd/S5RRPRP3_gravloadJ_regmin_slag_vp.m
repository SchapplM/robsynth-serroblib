% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = cos(pkin(8));
t23 = pkin(8) + qJ(4);
t18 = sin(t23);
t19 = cos(t23);
t32 = t19 * pkin(4) + t18 * qJ(5);
t44 = -t26 * pkin(3) - pkin(2) - t32;
t24 = qJ(1) + qJ(2);
t20 = sin(t24);
t21 = cos(t24);
t7 = g(1) * t20 - g(2) * t21;
t8 = g(1) * t21 + g(2) * t20;
t28 = sin(qJ(1));
t38 = t28 * pkin(1);
t27 = -pkin(7) - qJ(3);
t37 = t21 * t27;
t36 = t21 * pkin(2) + t20 * qJ(3);
t34 = t44 * t21;
t33 = -t20 * pkin(2) + t21 * qJ(3);
t30 = (-g(1) * t44 + g(2) * t27) * t20;
t29 = cos(qJ(1));
t22 = t29 * pkin(1);
t6 = t7 * t26;
t5 = t7 * sin(pkin(8));
t4 = t7 * t19;
t3 = t7 * t18;
t2 = g(3) * t18 + t8 * t19;
t1 = -g(3) * t19 + t8 * t18;
t9 = [0, g(1) * t28 - g(2) * t29, g(1) * t29 + g(2) * t28, 0, t7, t8, t6, -t5, -t8, -g(1) * (t33 - t38) - g(2) * (t22 + t36), 0, 0, 0, 0, 0, t4, -t3, t4, -t8, t3, -g(1) * (-t37 - t38) - g(2) * (t22 - t34) + t30; 0, 0, 0, 0, t7, t8, t6, -t5, -t8, -g(1) * t33 - g(2) * t36, 0, 0, 0, 0, 0, t4, -t3, t4, -t8, t3, g(1) * t37 + g(2) * t34 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t32 + t8 * (pkin(4) * t18 - qJ(5) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
