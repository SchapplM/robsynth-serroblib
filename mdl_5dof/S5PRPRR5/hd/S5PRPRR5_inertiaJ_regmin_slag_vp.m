% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t23 = sin(pkin(9));
t24 = cos(pkin(9));
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t14 = t26 * t23 - t29 * t24;
t19 = -t24 * pkin(3) - pkin(2);
t40 = 0.2e1 * t14 * pkin(4) + 0.2e1 * t19;
t39 = 0.2e1 * t19;
t25 = sin(qJ(5));
t38 = t25 * pkin(4);
t28 = cos(qJ(5));
t37 = t28 * pkin(4);
t36 = pkin(6) + qJ(3);
t35 = t23 ^ 2 + t24 ^ 2;
t16 = t36 * t23;
t17 = t36 * t24;
t34 = -t29 * t16 - t26 * t17;
t33 = t35 * qJ(3);
t32 = t26 * t16 - t29 * t17;
t15 = t29 * t23 + t26 * t24;
t30 = cos(qJ(2));
t27 = sin(qJ(2));
t11 = t14 * t27;
t10 = t15 * t27;
t8 = -t25 * t14 + t28 * t15;
t7 = t28 * t14 + t25 * t15;
t6 = -t14 * pkin(7) - t32;
t5 = -t15 * pkin(7) + t34;
t4 = t25 * t10 + t28 * t11;
t3 = -t28 * t10 + t25 * t11;
t2 = -t25 * t5 - t28 * t6;
t1 = -t25 * t6 + t28 * t5;
t9 = [1, 0, 0, 0, 0, 0, 0, t35 * t27 ^ 2 + t30 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t30, -t27, t30 * t24, -t30 * t23, t35 * t27, t30 * pkin(2) + t27 * t33, 0, 0, 0, 0, 0, -t30 * t14, -t30 * t15, 0, 0, 0, 0, 0, -t30 * t7, -t30 * t8; 0, 1, 0, 0, 0.2e1 * pkin(2) * t24, -0.2e1 * pkin(2) * t23, 0.2e1 * t33, t35 * qJ(3) ^ 2 + pkin(2) ^ 2, t15 ^ 2, -0.2e1 * t15 * t14, 0, 0, 0, t14 * t39, t15 * t39, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t40, t8 * t40; 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t24, t23, 0, -pkin(2), 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t11, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t34, t32, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
