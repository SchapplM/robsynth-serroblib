% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t26 = sin(pkin(9));
t28 = cos(pkin(9));
t36 = t26 ^ 2 + t28 ^ 2;
t37 = t36 * qJ(4);
t39 = t28 * pkin(4);
t29 = cos(pkin(8));
t19 = t29 * pkin(1) + pkin(2);
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t27 = sin(pkin(8));
t40 = pkin(1) * t27;
t10 = t33 * t19 - t31 * t40;
t9 = -pkin(3) - t10;
t4 = t9 - t39;
t44 = 0.2e1 * t4;
t20 = -pkin(3) - t39;
t43 = 0.2e1 * t20;
t42 = pkin(3) - t9;
t11 = -t31 * t19 - t33 * t40;
t8 = qJ(4) - t11;
t41 = t36 * t8;
t38 = t20 + t4;
t32 = cos(qJ(5));
t30 = sin(qJ(5));
t23 = t28 * pkin(7);
t16 = t28 * qJ(4) + t23;
t15 = (-pkin(7) - qJ(4)) * t26;
t14 = t32 * t26 + t30 * t28;
t13 = t30 * t26 - t32 * t28;
t12 = t14 ^ 2;
t3 = t28 * t8 + t23;
t2 = (-pkin(7) - t8) * t26;
t1 = -0.2e1 * t14 * t13;
t5 = [1, 0, 0, (t27 ^ 2 + t29 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t10, 0.2e1 * t11, -0.2e1 * t9 * t28, 0.2e1 * t9 * t26, 0.2e1 * t41, t36 * t8 ^ 2 + t9 ^ 2, t12, t1, 0, 0, 0, t13 * t44, t14 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t10, t11, t42 * t28, -t42 * t26, t37 + t41, -t9 * pkin(3) + t8 * t37, t12, t1, 0, 0, 0, t38 * t13, t38 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t28, -0.2e1 * pkin(3) * t26, 0.2e1 * t37, t36 * qJ(4) ^ 2 + pkin(3) ^ 2, t12, t1, 0, 0, 0, t13 * t43, t14 * t43; 0, 0, 0, 0, 0, 0, 0, -t28, t26, 0, t9, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t28, t26, 0, -pkin(3), 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t32 * t2 - t30 * t3, -t30 * t2 - t32 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t32 * t15 - t30 * t16, -t30 * t15 - t32 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
