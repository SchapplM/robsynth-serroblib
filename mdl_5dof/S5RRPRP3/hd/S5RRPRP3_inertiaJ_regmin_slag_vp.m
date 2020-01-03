% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t41 = t33 ^ 2 + t34 ^ 2;
t42 = t41 * qJ(3);
t35 = sin(qJ(4));
t45 = t35 * t33;
t46 = cos(qJ(4));
t17 = -t34 * t46 + t45;
t57 = 0.2e1 * t17;
t40 = t46 * t33;
t18 = t35 * t34 + t40;
t56 = -0.2e1 * t18;
t26 = -t34 * pkin(3) - pkin(2);
t49 = cos(qJ(2)) * pkin(1);
t19 = t26 - t49;
t55 = 0.2e1 * t19;
t54 = 0.2e1 * t26;
t50 = sin(qJ(2)) * pkin(1);
t25 = qJ(3) + t50;
t30 = t34 * pkin(7);
t14 = t34 * t25 + t30;
t47 = -pkin(7) - t25;
t6 = t35 * t14 - t40 * t47;
t7 = t14 * t46 + t45 * t47;
t53 = -t7 * t17 + t6 * t18;
t20 = t34 * qJ(3) + t30;
t39 = (-pkin(7) - qJ(3)) * t33;
t11 = t35 * t20 - t39 * t46;
t12 = t20 * t46 + t35 * t39;
t52 = t11 * t18 - t12 * t17;
t8 = t17 * pkin(4) - t18 * qJ(5) + t26;
t5 = t8 - t49;
t51 = t5 + t8;
t27 = -pkin(2) - t49;
t48 = pkin(2) - t27;
t44 = t19 + t26;
t43 = t41 * t25;
t15 = t18 ^ 2;
t10 = t17 * t56;
t9 = -pkin(4) * t18 - t17 * qJ(5);
t1 = [1, 0, 0, 1, 0.2e1 * t49, -0.2e1 * t50, -0.2e1 * t27 * t34, 0.2e1 * t27 * t33, 0.2e1 * t43, t25 ^ 2 * t41 + t27 ^ 2, t15, t10, 0, 0, 0, t17 * t55, t18 * t55, t5 * t57, 0.2e1 * t53, t5 * t56, t5 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 1, t49, -t50, t48 * t34, -t48 * t33, t42 + t43, -t27 * pkin(2) + t25 * t42, t15, t10, 0, 0, 0, t44 * t17, t44 * t18, t51 * t17, t52 + t53, -t51 * t18, t6 * t11 + t7 * t12 + t5 * t8; 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t34, -0.2e1 * pkin(2) * t33, 0.2e1 * t42, qJ(3) ^ 2 * t41 + pkin(2) ^ 2, t15, t10, 0, 0, 0, t17 * t54, t18 * t54, t8 * t57, 0.2e1 * t52, t8 * t56, t11 ^ 2 + t12 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, -t34, t33, 0, t27, 0, 0, 0, 0, 0, t17, t18, t17, 0, -t18, t5; 0, 0, 0, 0, 0, 0, -t34, t33, 0, -pkin(2), 0, 0, 0, 0, 0, t17, t18, t17, 0, -t18, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, -t6, -t7, -t6, t9, t7, -t6 * pkin(4) + t7 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, -t11, -t12, -t11, t9, t12, -t11 * pkin(4) + t12 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
