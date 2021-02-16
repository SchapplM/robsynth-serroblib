% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:55
% EndTime: 2021-01-15 23:23:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (195->52), mult. (255->86), div. (0->0), fcn. (275->10), ass. (0->46)
t30 = sin(qJ(3));
t33 = cos(qJ(3));
t35 = cos(qJ(1));
t41 = t35 * t33;
t32 = sin(qJ(1));
t34 = cos(qJ(2));
t47 = t32 * t34;
t15 = t30 * t47 + t41;
t42 = t35 * t30;
t17 = t32 * t33 - t34 * t42;
t31 = sin(qJ(2));
t49 = g(3) * t31;
t52 = -g(1) * t17 + g(2) * t15 + t30 * t49;
t38 = g(1) * t35 + g(2) * t32;
t11 = -g(3) * t34 + t38 * t31;
t28 = qJ(3) + pkin(9);
t27 = qJ(5) + t28;
t21 = sin(t27);
t46 = t35 * t21;
t22 = cos(t27);
t45 = t35 * t22;
t25 = sin(t28);
t44 = t35 * t25;
t26 = cos(t28);
t43 = t35 * t26;
t24 = t33 * pkin(3) + pkin(2);
t29 = qJ(4) + pkin(7);
t39 = t24 * t34 + t29 * t31;
t37 = g(1) * t32 - g(2) * t35;
t23 = t30 * pkin(3) + pkin(6);
t19 = t37 * t31;
t18 = t32 * t30 + t34 * t41;
t16 = -t33 * t47 + t42;
t13 = pkin(1) + t39;
t12 = t38 * t34 + t49;
t10 = t32 * t25 + t34 * t43;
t9 = t32 * t26 - t34 * t44;
t8 = -t26 * t47 + t44;
t7 = t25 * t47 + t43;
t6 = t32 * t21 + t34 * t45;
t5 = t32 * t22 - t34 * t46;
t4 = -t22 * t47 + t46;
t3 = t21 * t47 + t45;
t2 = g(1) * t6 - g(2) * t4 + t22 * t49;
t1 = -g(1) * t5 + g(2) * t3 + t21 * t49;
t14 = [0, t37, t38, 0, 0, 0, 0, 0, t37 * t34, -t19, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t19, -g(1) * (-t13 * t32 + t23 * t35) - g(2) * (t13 * t35 + t23 * t32), 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t11 * t33, -t11 * t30, t11 * t26, -t11 * t25, -t12, -g(3) * t39 - t38 * (-t31 * t24 + t29 * t34), 0, 0, 0, 0, 0, t11 * t22, -t11 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, g(1) * t18 - g(2) * t16 + t33 * t49, -g(1) * t9 + g(2) * t7 + t25 * t49, g(1) * t10 - g(2) * t8 + t26 * t49, 0, t52 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t14;
