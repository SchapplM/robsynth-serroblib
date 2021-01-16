% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:49
% EndTime: 2021-01-15 22:24:51
% DurationCPUTime: 0.34s
% Computational Cost: add. (518->57), mult. (975->102), div. (0->0), fcn. (1080->6), ass. (0->39)
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t36 = cos(qJ(2));
t48 = cos(qJ(3));
t19 = t34 * t35 - t48 * t36;
t20 = t34 * t36 + t48 * t35;
t29 = -t36 * pkin(2) - pkin(1);
t13 = t19 * pkin(3) + t29;
t54 = 0.2e1 * t13;
t53 = 0.2e1 * t29;
t52 = 0.2e1 * t36;
t51 = pkin(6) + pkin(7);
t33 = cos(pkin(8));
t50 = t33 * pkin(3);
t49 = t34 * pkin(2);
t31 = t48 * pkin(2);
t28 = t31 + pkin(3);
t32 = sin(pkin(8));
t17 = t32 * t28 + t33 * t49;
t11 = t20 * t51;
t38 = -t20 * qJ(4) - t11;
t12 = t19 * t51;
t8 = -t19 * qJ(4) - t12;
t4 = t32 * t8 - t33 * t38;
t6 = t32 * t38 + t33 * t8;
t45 = t4 ^ 2 + t6 ^ 2;
t30 = t32 * pkin(3);
t44 = -t30 - t17;
t41 = -t33 * t28 + t32 * t49;
t10 = -t32 * t19 + t33 * t20;
t9 = t33 * t19 + t32 * t20;
t40 = 0.2e1 * t4 * t10 - 0.2e1 * t6 * t9;
t39 = -t41 + t50;
t25 = pkin(4) + t50;
t24 = t30 + qJ(5);
t15 = -pkin(4) + t41;
t14 = qJ(5) + t17;
t2 = t9 * pkin(4) - t10 * qJ(5) + t13;
t1 = [1, 0, 0, t35 ^ 2, t35 * t52, 0, 0, 0, pkin(1) * t52, -0.2e1 * pkin(1) * t35, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t53, t20 * t53, t9 * t54, t10 * t54, t40, t13 ^ 2 + t45, 0.2e1 * t2 * t9, t40, -0.2e1 * t2 * t10, t2 ^ 2 + t45; 0, 0, 0, 0, 0, t35, t36, 0, -t35 * pkin(6), -t36 * pkin(6), 0, 0, t20, -t19, 0, -t11, t12, -t4, -t6, t10 * t41 - t17 * t9, t6 * t17 + t4 * t41, -t4, t15 * t10 - t14 * t9, t6, t6 * t14 + t4 * t15; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t49, -0.2e1 * t41, -0.2e1 * t17, 0, t17 ^ 2 + t41 ^ 2, -0.2e1 * t15, 0, 0.2e1 * t14, t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, -t11, t12, -t4, -t6, (-t10 * t33 - t32 * t9) * pkin(3), (t32 * t6 - t33 * t4) * pkin(3), -t4, -t25 * t10 - t24 * t9, t6, t6 * t24 - t4 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t49, t39, t44, 0, (t17 * t32 - t33 * t41) * pkin(3), 0.2e1 * pkin(4) + t39, 0, 0.2e1 * qJ(5) - t44, t14 * t24 - t15 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t50, -0.2e1 * t30, 0, (t32 ^ 2 + t33 ^ 2) * pkin(3) ^ 2, 0.2e1 * t25, 0, 0.2e1 * t24, t24 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t13, t9, 0, -t10, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
