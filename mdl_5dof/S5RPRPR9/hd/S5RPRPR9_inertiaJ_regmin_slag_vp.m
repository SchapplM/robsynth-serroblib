% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t38 = 2 * qJ(4);
t23 = -pkin(3) - pkin(7);
t22 = cos(qJ(3));
t37 = t22 * pkin(3);
t19 = sin(qJ(5));
t20 = sin(qJ(3));
t36 = t19 * t20;
t35 = t19 * t22;
t34 = t20 * t22;
t21 = cos(qJ(5));
t33 = t21 * t19;
t32 = t21 * t22;
t14 = t20 ^ 2;
t16 = t22 ^ 2;
t31 = t14 + t16;
t30 = t20 * qJ(4);
t29 = t22 * qJ(4);
t28 = -0.2e1 * t34;
t18 = cos(pkin(8));
t11 = -t18 * pkin(1) - pkin(2);
t27 = -t20 * pkin(3) + t29;
t26 = t20 * t23 + t29;
t25 = t11 - t30;
t17 = sin(pkin(8));
t15 = t21 ^ 2;
t13 = t19 ^ 2;
t12 = t21 * t20;
t10 = t17 * pkin(1) + pkin(6);
t8 = t22 * t10;
t7 = t20 * t10;
t6 = t22 * pkin(4) + t8;
t5 = t20 * pkin(4) + t7;
t4 = t25 - t37;
t3 = t23 * t22 + t25;
t2 = t19 * t5 + t21 * t3;
t1 = -t19 * t3 + t21 * t5;
t9 = [1, 0, 0, (t17 ^ 2 + t18 ^ 2) * pkin(1) ^ 2, t14, 0.2e1 * t34, 0, 0, 0, -0.2e1 * t11 * t22, 0.2e1 * t11 * t20, 0.2e1 * t31 * t10, 0.2e1 * t4 * t22, -0.2e1 * t4 * t20, t31 * t10 ^ 2 + t4 ^ 2, t13 * t16, 0.2e1 * t16 * t33, t19 * t28, t21 * t28, t14, 0.2e1 * t1 * t20 + 0.2e1 * t6 * t32, -0.2e1 * t2 * t20 - 0.2e1 * t6 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t20, t22, 0, -t7, -t8, t27, t7, t8, t27 * t10, -t19 * t32, (t13 - t15) * t22, t12, -t36, 0, t6 * t19 + t26 * t21, -t26 * t19 + t6 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, -t22, t20, t30 + t37, 0, 0, 0, 0, 0, t36, t12; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t38, pkin(3) ^ 2 + (qJ(4) ^ 2), t15, -0.2e1 * t33, 0, 0, 0, t19 * t38, t21 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, t7, 0, 0, 0, 0, 0, t12, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t32, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, t21 * t23, -t19 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
