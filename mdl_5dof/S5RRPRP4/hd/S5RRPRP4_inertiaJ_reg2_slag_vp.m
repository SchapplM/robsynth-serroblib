% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP4
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:55
% EndTime: 2019-12-31 19:52:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (183->47), mult. (284->62), div. (0->0), fcn. (201->4), ass. (0->42)
t29 = sin(qJ(4));
t26 = t29 ^ 2;
t31 = cos(qJ(4));
t27 = t31 ^ 2;
t14 = t26 + t27;
t33 = -pkin(2) - pkin(7);
t8 = t14 * t33;
t30 = sin(qJ(2));
t47 = t30 * pkin(1);
t19 = qJ(3) + t47;
t53 = t19 ^ 2;
t52 = 0.2e1 * t19;
t51 = 0.2e1 * t29;
t50 = -0.2e1 * t31;
t34 = 0.2e1 * qJ(3);
t9 = t29 * pkin(4) - t31 * qJ(5) + qJ(3);
t4 = t9 + t47;
t49 = t4 + t9;
t32 = cos(qJ(2));
t46 = t32 * pkin(1);
t21 = -pkin(2) - t46;
t16 = -pkin(7) + t21;
t48 = t16 * t8;
t44 = t29 * t16;
t43 = t29 * t33;
t42 = t31 * t29;
t41 = t14 * t16 ^ 2;
t40 = t14 * t33 ^ 2;
t39 = t19 * qJ(3);
t38 = qJ(3) + t19;
t13 = t31 * pkin(4) + t29 * qJ(5);
t36 = -0.2e1 * pkin(2);
t35 = qJ(3) ^ 2;
t24 = t31 * t33;
t18 = -0.2e1 * t42;
t17 = 0.2e1 * t42;
t12 = t31 * t16;
t5 = 0.2e1 * t8;
t3 = t14 * t16;
t2 = 0.2e1 * t3;
t1 = (-t16 - t33) * t14;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t47, 0, (t30 ^ 2 + t32 ^ 2) * pkin(1) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t21, t52, t21 ^ 2 + t53, t27, t18, 0, t26, 0, 0, t19 * t51, t31 * t52, -t2, t41 + t53, t27, 0, t17, 0, 0, t26, t4 * t51, -t2, t4 * t50, t4 ^ 2 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t47, 0, 0, 1, 0, 0, 0, 0, 0, 0, t36 - t46, t34 + t47, -t21 * pkin(2) + t39, t27, t18, 0, t26, 0, 0, t38 * t29, t38 * t31, t1, t39 + t48, t27, 0, t17, 0, 0, t26, t49 * t29, t1, -t49 * t31, t4 * t9 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t36, t34, pkin(2) ^ 2 + t35, t27, t18, 0, t26, 0, 0, t29 * t34, t31 * t34, -t5, t35 + t40, t27, 0, t17, 0, 0, t26, t9 * t51, -t5, t9 * t50, t9 ^ 2 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t3, 0, 0, 0, 0, 0, 0, 0, -t14, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t14, t8, 0, 0, 0, 0, 0, 0, 0, -t14, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, 0, t12, -t44, 0, 0, 0, t31, 0, 0, t29, 0, t12, -t13, t44, t13 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, 0, t24, -t43, 0, 0, 0, t31, 0, 0, t29, 0, t24, -t13, t43, t13 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t29, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
