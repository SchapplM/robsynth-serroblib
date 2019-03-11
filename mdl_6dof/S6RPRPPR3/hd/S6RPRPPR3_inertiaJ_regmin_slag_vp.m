% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t30 = sin(qJ(3));
t48 = -0.2e1 * t30;
t47 = 0.2e1 * t30;
t32 = cos(qJ(3));
t46 = -0.2e1 * t32;
t33 = -pkin(3) - pkin(4);
t27 = cos(pkin(9));
t13 = -t27 * pkin(1) - pkin(2);
t29 = sin(qJ(6));
t45 = t29 * t30;
t31 = cos(qJ(6));
t44 = t29 * t31;
t14 = t29 * t32;
t43 = t31 * t32;
t26 = sin(pkin(9));
t12 = t26 * pkin(1) + pkin(7);
t42 = t32 * t12;
t16 = t30 * qJ(4);
t41 = t32 * pkin(3) + t16;
t23 = t30 ^ 2;
t25 = t32 ^ 2;
t10 = t23 + t25;
t40 = t32 * qJ(4);
t39 = t32 * t47;
t5 = -t41 + t13;
t4 = t32 * pkin(4) - t5;
t38 = -t30 * pkin(3) + t40;
t21 = -pkin(8) + t33;
t28 = qJ(4) + pkin(5);
t37 = -t21 * t30 - t28 * t32;
t35 = qJ(4) ^ 2;
t34 = 0.2e1 * qJ(4);
t24 = t31 ^ 2;
t22 = t29 ^ 2;
t15 = t31 * t30;
t9 = t30 * t12;
t7 = -t32 * qJ(5) + t42;
t6 = -t30 * qJ(5) + t9;
t3 = t30 * pkin(5) + t32 * pkin(8) + t4;
t2 = t29 * t3 + t31 * t6;
t1 = -t29 * t6 + t31 * t3;
t8 = [1, 0, 0 (t26 ^ 2 + t27 ^ 2) * pkin(1) ^ 2, t23, t39, 0, 0, 0, t13 * t46, t13 * t47, t5 * t46, 0.2e1 * t10 * t12, t5 * t48, t10 * t12 ^ 2 + t5 ^ 2, t4 * t47, t4 * t46, -0.2e1 * t6 * t30 - 0.2e1 * t7 * t32, t4 ^ 2 + t6 ^ 2 + t7 ^ 2, t24 * t25, -0.2e1 * t25 * t44, t43 * t48, t29 * t39, t23, 0.2e1 * t1 * t30 - 0.2e1 * t7 * t14, -0.2e1 * t2 * t30 - 0.2e1 * t7 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t30 - t6 * t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t30, t32, 0, -t9, -t42, -t9, t38, t42, t38 * t12, t7, t6, -t33 * t30 - t40, t7 * qJ(4) + t6 * t33, t29 * t43 (-t22 + t24) * t32, -t45, -t15, 0, t29 * t37 + t7 * t31, -t7 * t29 + t31 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, t32, 0, t30, t41, t30, -t32, 0, -t32 * t33 + t16, 0, 0, 0, 0, 0, t15, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t34, pkin(3) ^ 2 + t35, t34, 0.2e1 * t33, 0, t33 ^ 2 + t35, t22, 0.2e1 * t44, 0, 0, 0, 0.2e1 * t28 * t31, -0.2e1 * t28 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t9, 0, 0, -t30, t6, 0, 0, 0, 0, 0, -t45, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 1, 0, t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t32, 0, t4, 0, 0, 0, 0, 0, t15, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t14, t30, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t31, 0, -t29 * t21, -t31 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
