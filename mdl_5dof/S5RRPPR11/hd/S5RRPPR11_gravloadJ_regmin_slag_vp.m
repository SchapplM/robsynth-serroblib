% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:00
% DurationCPUTime: 0.21s
% Computational Cost: add. (107->53), mult. (215->73), div. (0->0), fcn. (217->8), ass. (0->40)
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t35 = t27 * pkin(2) + t25 * qJ(3);
t30 = -pkin(1) - t35;
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t9 = g(1) * t28 + g(2) * t26;
t51 = t9 * t25;
t6 = g(3) * t25 + t9 * t27;
t50 = pkin(2) * t25;
t49 = g(1) * t26;
t45 = g(3) * t27;
t22 = pkin(8) + qJ(5);
t14 = sin(t22);
t44 = t26 * t14;
t15 = cos(t22);
t43 = t26 * t15;
t23 = sin(pkin(8));
t42 = t26 * t23;
t24 = cos(pkin(8));
t41 = t26 * t24;
t40 = t28 * t14;
t39 = t28 * t15;
t38 = t28 * t23;
t37 = t28 * t24;
t34 = qJ(3) * t27;
t33 = t27 * qJ(4);
t32 = t26 * pkin(6) - t30 * t28;
t31 = -g(2) * t28 + t49;
t19 = t28 * pkin(6);
t12 = t28 * t34;
t10 = t26 * t34;
t8 = t31 * t27;
t7 = t31 * t25;
t5 = -t45 + t51;
t4 = -t25 * t44 + t39;
t3 = t25 * t43 + t40;
t2 = t25 * t40 + t43;
t1 = t25 * t39 - t44;
t11 = [0, t31, t9, 0, 0, 0, 0, 0, t8, -t7, -t9, -t8, t7, -g(1) * t19 - g(2) * t32 - t30 * t49, -g(1) * (-t25 * t42 + t37) - g(2) * (t25 * t38 + t41), -g(1) * (-t25 * t41 - t38) - g(2) * (t25 * t37 - t42), t8, -g(1) * (t28 * pkin(3) + t19) - g(2) * (t28 * t33 + t32) + (-g(1) * (t30 - t33) - g(2) * pkin(3)) * t26, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, -t5, -t6, -g(1) * (-t28 * t50 + t12) - g(2) * (-t26 * t50 + t10) - g(3) * t35, -t6 * t23, -t6 * t24, t5, -g(1) * t12 - g(2) * t10 - g(3) * (t33 + t35) + (pkin(2) + qJ(4)) * t51, 0, 0, 0, 0, 0, -t6 * t14, -t6 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t15 * t45, g(1) * t2 - g(2) * t4 - t14 * t45;];
taug_reg = t11;
