% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:12
% DurationCPUTime: 0.27s
% Computational Cost: add. (350->46), mult. (667->79), div. (0->0), fcn. (779->6), ass. (0->33)
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t24 = sin(qJ(3));
t25 = cos(qJ(3));
t10 = t24 * t21 - t25 * t22;
t11 = t25 * t21 + t24 * t22;
t23 = sin(qJ(4));
t35 = cos(qJ(4));
t7 = -t23 * t10 + t35 * t11;
t38 = -0.2e1 * t7;
t15 = -t22 * pkin(2) - pkin(1);
t8 = t10 * pkin(3) + t15;
t37 = 0.2e1 * t8;
t36 = 0.2e1 * t15;
t34 = pkin(6) + qJ(2);
t33 = t21 ^ 2 + t22 ^ 2;
t32 = t35 * pkin(3);
t12 = t34 * t21;
t13 = t34 * t22;
t31 = -t25 * t12 - t24 * t13;
t30 = t24 * t12 - t25 * t13;
t29 = -t11 * pkin(7) + t31;
t28 = 2 * pkin(4);
t26 = 2 * qJ(5);
t18 = t23 * pkin(3);
t16 = t32 + pkin(4);
t14 = t18 + qJ(5);
t6 = t35 * t10 + t23 * t11;
t5 = -t10 * pkin(7) - t30;
t3 = t23 * t29 + t35 * t5;
t2 = t23 * t5 - t35 * t29;
t1 = t6 * pkin(4) - t7 * qJ(5) + t8;
t4 = [1, 0, 0, 0.2e1 * pkin(1) * t22, -0.2e1 * pkin(1) * t21, 0.2e1 * t33 * qJ(2), t33 * qJ(2) ^ 2 + pkin(1) ^ 2, t11 ^ 2, -0.2e1 * t11 * t10, 0, 0, 0, t10 * t36, t11 * t36, t7 ^ 2, t6 * t38, 0, 0, 0, t6 * t37, t7 * t37, 0.2e1 * t1 * t6, 0.2e1 * t2 * t7 - 0.2e1 * t3 * t6, t1 * t38, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, -t22, t21, 0, -pkin(1), 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, t6, t7, t6, 0, -t7, t1; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, t31, t30, 0, 0, t7, -t6, 0, -t2, -t3, -t2, -t14 * t6 - t16 * t7, t3, t3 * t14 - t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t18, 0.2e1 * t16, 0, 0.2e1 * t14, t14 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, -t2, -t3, -t2, -pkin(4) * t7 - t6 * qJ(5), t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t18, t28 + t32, 0, t26 + t18, t16 * pkin(4) + t14 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, 0, t26, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
