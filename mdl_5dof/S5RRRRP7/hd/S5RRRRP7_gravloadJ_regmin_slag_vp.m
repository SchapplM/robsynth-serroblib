% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:44
% DurationCPUTime: 0.22s
% Computational Cost: add. (233->54), mult. (328->78), div. (0->0), fcn. (343->8), ass. (0->43)
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t55 = pkin(4) * t29 + qJ(5) * t26;
t25 = qJ(2) + qJ(3);
t22 = sin(t25);
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t39 = g(1) * t31 + g(2) * t28;
t35 = t39 * t22;
t23 = cos(t25);
t54 = t23 * pkin(3) + t22 * pkin(8);
t27 = sin(qJ(2));
t53 = pkin(2) * t27;
t51 = pkin(8) * t23;
t48 = g(3) * t22;
t47 = g(3) * t26;
t46 = t28 * t26;
t45 = t28 * t29;
t44 = t31 * t26;
t43 = t31 * t29;
t41 = t55 * t23 + t54;
t10 = t23 * t44 - t45;
t8 = t23 * t46 + t43;
t40 = g(1) * t8 - g(2) * t10;
t38 = g(1) * t28 - g(2) * t31;
t30 = cos(qJ(2));
t24 = t30 * pkin(2);
t37 = t24 + pkin(1) + t54;
t1 = g(1) * t10 + g(2) * t8 + t22 * t47;
t11 = t23 * t43 + t46;
t9 = t23 * t45 - t44;
t34 = g(1) * t11 + g(2) * t9 + t29 * t48;
t5 = -g(3) * t23 + t35;
t33 = (pkin(3) + t55) * t35;
t32 = -pkin(7) - pkin(6);
t15 = t31 * t51;
t13 = t28 * t51;
t7 = t38 * t22;
t6 = t23 * t39 + t48;
t4 = t5 * t29;
t3 = -t23 * t47 + t26 * t35;
t2 = g(1) * t9 - g(2) * t11;
t12 = [0, t38, t39, 0, 0, 0, 0, 0, t38 * t30, -t38 * t27, 0, 0, 0, 0, 0, t38 * t23, -t7, 0, 0, 0, 0, 0, t2, -t40, t2, t7, t40, -g(1) * (-t9 * pkin(4) - t8 * qJ(5)) - g(2) * (t11 * pkin(4) + t10 * qJ(5)) + (g(1) * t32 - g(2) * t37) * t31 + (g(1) * t37 + g(2) * t32) * t28; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t30 + t27 * t39, g(3) * t27 + t30 * t39, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-t31 * t53 + t15) - g(2) * (-t28 * t53 + t13) - g(3) * (t24 + t41) + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t15 - g(2) * t13 - g(3) * t41 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t34, t1, 0, -t34, -g(1) * (-t10 * pkin(4) + t11 * qJ(5)) - g(2) * (-t8 * pkin(4) + t9 * qJ(5)) - (-pkin(4) * t26 + qJ(5) * t29) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
