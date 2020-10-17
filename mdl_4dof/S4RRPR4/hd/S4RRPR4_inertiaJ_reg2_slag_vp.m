% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.30s
% Computational Cost: add. (183->40), mult. (370->77), div. (0->0), fcn. (364->6), ass. (0->40)
t33 = sin(pkin(7));
t31 = t33 ^ 2;
t34 = cos(pkin(7));
t32 = t34 ^ 2;
t40 = t31 + t32;
t41 = t40 * qJ(3);
t26 = -t34 * pkin(3) - pkin(2);
t37 = cos(qJ(2));
t46 = t37 * pkin(1);
t17 = t26 - t46;
t52 = 0.2e1 * t17;
t51 = 0.2e1 * t26;
t50 = 0.2e1 * t34;
t35 = sin(qJ(4));
t44 = cos(qJ(4));
t14 = t35 * t33 - t44 * t34;
t16 = t44 * t33 + t35 * t34;
t36 = sin(qJ(2));
t47 = t36 * pkin(1);
t25 = qJ(3) + t47;
t10 = (-pkin(6) - t25) * t33;
t30 = t34 * pkin(6);
t11 = t34 * t25 + t30;
t5 = t44 * t10 - t35 * t11;
t6 = t35 * t10 + t44 * t11;
t49 = -t6 * t14 - t5 * t16;
t18 = (-pkin(6) - qJ(3)) * t33;
t19 = t34 * qJ(3) + t30;
t8 = t44 * t18 - t35 * t19;
t9 = t35 * t18 + t44 * t19;
t48 = -t9 * t14 - t8 * t16;
t27 = -pkin(2) - t46;
t45 = pkin(2) - t27;
t43 = t17 + t26;
t42 = t40 * t25;
t22 = t33 * t50;
t13 = t16 ^ 2;
t12 = t14 ^ 2;
t7 = -0.2e1 * t16 * t14;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t47, 0, (t36 ^ 2 + t37 ^ 2) * pkin(1) ^ 2, t31, t22, 0, t32, 0, 0, -0.2e1 * t27 * t34, 0.2e1 * t27 * t33, 0.2e1 * t42, t40 * t25 ^ 2 + t27 ^ 2, t13, t7, 0, t12, 0, 0, t14 * t52, t16 * t52, 0.2e1 * t49, t17 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t47, 0, 0, t31, t22, 0, t32, 0, 0, t45 * t34, -t45 * t33, t41 + t42, -t27 * pkin(2) + t25 * t41, t13, t7, 0, t12, 0, 0, t43 * t14, t43 * t16, t48 + t49, t17 * t26 + t5 * t8 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, t22, 0, t32, 0, 0, pkin(2) * t50, -0.2e1 * pkin(2) * t33, 0.2e1 * t41, t40 * qJ(3) ^ 2 + pkin(2) ^ 2, t13, t7, 0, t12, 0, 0, t14 * t51, t16 * t51, 0.2e1 * t48, t26 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, t27, 0, 0, 0, 0, 0, 0, t14, t16, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t14, t16, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
