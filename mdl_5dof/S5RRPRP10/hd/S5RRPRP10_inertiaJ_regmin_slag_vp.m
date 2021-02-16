% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:44
% EndTime: 2021-01-15 21:02:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (232->65), mult. (412->119), div. (0->0), fcn. (384->4), ass. (0->48)
t26 = sin(qJ(4));
t16 = t26 * pkin(4) + qJ(3);
t51 = 0.2e1 * t16;
t27 = sin(qJ(2));
t50 = -0.2e1 * t27;
t29 = cos(qJ(2));
t49 = 0.2e1 * t29;
t48 = 0.2e1 * qJ(3);
t30 = -pkin(2) - pkin(7);
t47 = t27 * pkin(4);
t28 = cos(qJ(4));
t46 = t28 * pkin(4);
t19 = t27 * pkin(6);
t12 = t27 * pkin(3) + t19;
t45 = t26 * t12;
t44 = t26 * t27;
t43 = t26 * t29;
t42 = t27 * t29;
t41 = t28 * t26;
t40 = t28 * t29;
t20 = t29 * pkin(6);
t13 = t29 * pkin(3) + t20;
t23 = t27 ^ 2;
t25 = t29 ^ 2;
t39 = t23 + t25;
t38 = qJ(5) * t29;
t37 = t29 * qJ(3);
t36 = -0.2e1 * t42;
t35 = -t27 * qJ(3) - pkin(1);
t7 = t30 * t29 + t35;
t3 = t28 * t12 - t26 * t7;
t34 = t26 * t38 + t3;
t33 = -t27 * pkin(2) + t37;
t32 = t27 * t30 + t37;
t24 = t28 ^ 2;
t22 = t26 ^ 2;
t18 = t28 * t30;
t17 = t28 * t27;
t14 = t22 + t24;
t11 = -t29 * pkin(2) + t35;
t10 = -t28 * qJ(5) + t18;
t9 = (-qJ(5) + t30) * t26;
t6 = pkin(4) * t40 + t13;
t5 = t10 * t28 + t9 * t26;
t4 = t28 * t7 + t45;
t2 = t45 + (t7 - t38) * t28;
t1 = t34 + t47;
t8 = [1, 0, 0, t23, 0.2e1 * t42, 0, 0, 0, pkin(1) * t49, pkin(1) * t50, 0.2e1 * t39 * pkin(6), t11 * t49, t11 * t50, t39 * pkin(6) ^ 2 + t11 ^ 2, t22 * t25, 0.2e1 * t25 * t41, t26 * t36, t28 * t36, t23, 0.2e1 * t13 * t40 + 0.2e1 * t3 * t27, -0.2e1 * t13 * t43 - 0.2e1 * t4 * t27, 0.2e1 * t1 * t27 + 0.2e1 * t6 * t40, -0.2e1 * t2 * t27 - 0.2e1 * t6 * t43, (t1 * t26 - t2 * t28) * t49, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t27, t29, 0, -t19, -t20, t33, t19, t20, t33 * pkin(6), -t26 * t40, (t22 - t24) * t29, t17, -t44, 0, t13 * t26 + t32 * t28, t13 * t28 - t32 * t26, t10 * t27 + t16 * t40 + t6 * t26, -t16 * t43 - t9 * t27 + t6 * t28, (-t29 * t9 - t1) * t28 + (t10 * t29 - t2) * t26, t1 * t10 + t6 * t16 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t48, pkin(2) ^ 2 + qJ(3) ^ 2, t24, -0.2e1 * t41, 0, 0, 0, t26 * t48, t28 * t48, t26 * t51, t28 * t51, -0.2e1 * t5, t10 ^ 2 + t16 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, t19, 0, 0, 0, 0, 0, t17, -t44, t17, -t44, 0, t1 * t28 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t40, t27, t3, -t4, t34 + 0.2e1 * t47, -t2, pkin(4) * t43, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, t18, -t26 * t30, t10, -t9, -t46, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, t28, -t26, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t43, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
