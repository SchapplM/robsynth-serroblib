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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:41
% EndTime: 2021-01-15 15:52:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (189->45), mult. (420->92), div. (0->0), fcn. (509->8), ass. (0->35)
t25 = sin(pkin(9));
t26 = cos(pkin(9));
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t17 = t25 * t28 - t26 * t31;
t24 = -t31 * pkin(3) - pkin(2);
t39 = 0.2e1 * t17 * pkin(4) + 0.2e1 * t24;
t38 = 0.2e1 * t24;
t37 = 0.2e1 * t31;
t36 = t25 * pkin(3);
t35 = t26 * pkin(3);
t34 = qJ(4) + pkin(6);
t20 = t34 * t28;
t21 = t34 * t31;
t9 = -t26 * t20 - t25 * t21;
t10 = -t25 * t20 + t26 * t21;
t18 = t25 * t31 + t26 * t28;
t32 = cos(qJ(2));
t30 = cos(qJ(5));
t29 = sin(qJ(2));
t27 = sin(qJ(5));
t23 = pkin(4) + t35;
t15 = -t27 * t23 - t30 * t36;
t14 = t30 * t23 - t27 * t36;
t13 = t17 * t29;
t12 = t18 * t29;
t8 = -t27 * t17 + t30 * t18;
t7 = t30 * t17 + t27 * t18;
t6 = -t17 * pkin(7) + t10;
t5 = -t18 * pkin(7) + t9;
t4 = t27 * t12 + t30 * t13;
t3 = -t30 * t12 + t27 * t13;
t2 = -t27 * t5 - t30 * t6;
t1 = -t27 * t6 + t30 * t5;
t11 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 ^ 2 + t13 ^ 2 + t32 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t32, -t29, 0, 0, 0, 0, 0, t32 * t31, -t32 * t28, -t32 * t17, -t32 * t18, t12 * t18 + t13 * t17, -t13 * t10 - t12 * t9 - t32 * t24, 0, 0, 0, 0, 0, -t32 * t7, -t32 * t8; 0, 1, 0, 0, t28 ^ 2, t28 * t37, 0, 0, 0, pkin(2) * t37, -0.2e1 * pkin(2) * t28, t17 * t38, t18 * t38, -0.2e1 * t10 * t17 - 0.2e1 * t9 * t18, t10 ^ 2 + t24 ^ 2 + t9 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t39, t8 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t29, -t31 * t29, -t12, t13, 0, (-t12 * t26 - t13 * t25) * pkin(3), 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(6), -t31 * pkin(6), t9, -t10, (-t17 * t25 - t18 * t26) * pkin(3), (t10 * t25 + t26 * t9) * pkin(3), 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t35, -0.2e1 * t36, 0, (t25 ^ 2 + t26 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t14, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, t24, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t11;
