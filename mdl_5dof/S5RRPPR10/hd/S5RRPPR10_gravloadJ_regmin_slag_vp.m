% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:40
% EndTime: 2019-12-31 19:44:41
% DurationCPUTime: 0.24s
% Computational Cost: add. (122->53), mult. (330->84), div. (0->0), fcn. (366->8), ass. (0->45)
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t46 = t31 * pkin(2) + t28 * qJ(3);
t57 = -pkin(1) - t46;
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t42 = g(1) * t32 + g(2) * t29;
t56 = t42 * t28;
t6 = -g(3) * t31 + t56;
t55 = pkin(2) * t28;
t54 = g(1) * t29;
t51 = g(3) * t28;
t49 = t29 * t31;
t25 = sin(pkin(8));
t48 = t32 * t25;
t26 = cos(pkin(8));
t47 = t32 * t26;
t45 = qJ(3) * t31;
t44 = t29 * pkin(6) - t57 * t32;
t10 = -t29 * t26 + t31 * t48;
t8 = t25 * t49 + t47;
t43 = g(1) * t8 - g(2) * t10;
t41 = -g(2) * t32 + t54;
t27 = sin(qJ(5));
t30 = cos(qJ(5));
t9 = t26 * t49 - t48;
t40 = t9 * t27 - t8 * t30;
t39 = t8 * t27 + t9 * t30;
t38 = pkin(3) * t26 + qJ(4) * t25;
t37 = t25 * t30 - t26 * t27;
t36 = t25 * t27 + t26 * t30;
t34 = g(3) * t37;
t33 = t57 * t54;
t22 = t32 * pkin(6);
t15 = t32 * t45;
t13 = t29 * t45;
t12 = t41 * t28;
t11 = t29 * t25 + t31 * t47;
t7 = t42 * t31 + t51;
t5 = t6 * t26;
t4 = t6 * t25;
t3 = g(1) * t9 - g(2) * t11;
t2 = t10 * t27 + t11 * t30;
t1 = t10 * t30 - t11 * t27;
t14 = [0, t41, t42, 0, 0, 0, 0, 0, t41 * t31, -t12, t3, -t43, t12, -g(1) * t22 - g(2) * t44 - t33, t3, t12, t43, -g(1) * (-t9 * pkin(3) - t8 * qJ(4) + t22) - g(2) * (t11 * pkin(3) + t10 * qJ(4) + t44) - t33, 0, 0, 0, 0, 0, g(1) * t39 - g(2) * t2, -g(1) * t40 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, t5, -t4, -t7, -g(1) * (-t32 * t55 + t15) - g(2) * (-t29 * t55 + t13) - g(3) * t46, t5, -t7, t4, -g(1) * t15 - g(2) * t13 - g(3) * (t38 * t31 + t46) + (pkin(2) + t38) * t56, 0, 0, 0, 0, 0, t6 * t36, -t31 * t34 + t37 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8 - t25 * t51, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t40 - t28 * t34, g(1) * t2 + g(2) * t39 + t36 * t51;];
taug_reg = t14;
