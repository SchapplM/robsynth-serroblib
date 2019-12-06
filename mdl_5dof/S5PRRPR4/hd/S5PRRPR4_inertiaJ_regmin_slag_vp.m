% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t26 = sin(pkin(9));
t27 = cos(pkin(9));
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t17 = -t26 * t29 + t27 * t32;
t25 = -t32 * pkin(3) - pkin(2);
t38 = -0.2e1 * t17 * pkin(4) + 0.2e1 * t25;
t37 = 0.2e1 * t32;
t36 = pkin(3) * t26;
t35 = -qJ(4) - pkin(6);
t22 = t35 * t29;
t23 = t35 * t32;
t10 = t26 * t22 - t27 * t23;
t9 = t27 * t22 + t26 * t23;
t18 = t26 * t32 + t27 * t29;
t33 = cos(qJ(2));
t31 = cos(qJ(5));
t30 = sin(qJ(2));
t28 = sin(qJ(5));
t24 = t27 * pkin(3) + pkin(4);
t15 = -t28 * t24 - t31 * t36;
t14 = t31 * t24 - t28 * t36;
t13 = t17 * t30;
t12 = t18 * t30;
t8 = t28 * t17 + t31 * t18;
t7 = -t31 * t17 + t28 * t18;
t6 = t17 * pkin(7) + t10;
t5 = -t18 * pkin(7) + t9;
t4 = t28 * t12 - t31 * t13;
t3 = -t31 * t12 - t28 * t13;
t2 = -t28 * t5 - t31 * t6;
t1 = -t28 * t6 + t31 * t5;
t11 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 ^ 2 + t13 ^ 2 + t33 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t33, -t30, 0, 0, 0, 0, 0, t33 * t32, -t33 * t29, t12 * t18 + t13 * t17, t13 * t10 - t12 * t9 - t33 * t25, 0, 0, 0, 0, 0, -t33 * t7, -t33 * t8; 0, 1, 0, 0, t29 ^ 2, t29 * t37, 0, 0, 0, pkin(2) * t37, -0.2e1 * pkin(2) * t29, 0.2e1 * t10 * t17 - 0.2e1 * t9 * t18, t10 ^ 2 + t25 ^ 2 + t9 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t38, t8 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29 * t30, -t32 * t30, 0, (-t12 * t27 + t13 * t26) * pkin(3), 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, t29, t32, 0, -t29 * pkin(6), -t32 * pkin(6), (t17 * t26 - t18 * t27) * pkin(3), (t10 * t26 + t27 * t9) * pkin(3), 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t26 ^ 2 + t27 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t14, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t11;
