% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:15:56
% EndTime: 2019-05-05 14:15:57
% DurationCPUTime: 0.43s
% Computational Cost: add. (374->71), mult. (606->132), div. (0->0), fcn. (709->8), ass. (0->43)
t38 = sin(pkin(10));
t40 = cos(pkin(10));
t43 = sin(qJ(4));
t56 = cos(qJ(4));
t21 = t38 * t56 + t40 * t43;
t18 = t21 ^ 2;
t42 = sin(qJ(6));
t47 = -t38 * t43 + t40 * t56;
t16 = t47 * t42;
t44 = cos(qJ(6));
t15 = t47 * t44;
t55 = t42 * t21;
t54 = t42 * t44;
t53 = t44 * t21;
t39 = sin(pkin(9));
t41 = cos(pkin(9));
t45 = -pkin(1) - pkin(2);
t24 = t39 * qJ(2) - t41 * t45;
t26 = t41 * qJ(2) + t39 * t45;
t52 = 0.2e1 * t56;
t23 = pkin(3) + t24;
t51 = -pkin(7) + t26;
t14 = t56 * pkin(4) + t23;
t28 = t38 * pkin(4) + pkin(8);
t29 = -t40 * pkin(4) - pkin(5);
t50 = -t21 * t29 - t28 * t47;
t49 = (qJ(5) - t51) * t43;
t48 = t56 * t51;
t37 = t44 ^ 2;
t36 = t42 ^ 2;
t35 = t41 ^ 2;
t17 = t47 ^ 2;
t13 = t47 * t39;
t11 = t21 * t39;
t10 = -t56 * qJ(5) + t48;
t8 = t44 * t13 - t42 * t41;
t7 = -t42 * t13 - t44 * t41;
t6 = pkin(5) * t47 + t21 * pkin(8) + t14;
t5 = t40 * t10 + t38 * t49;
t3 = t38 * t10 - t40 * t49;
t2 = t42 * t6 + t44 * t5;
t1 = -t42 * t5 + t44 * t6;
t4 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t24, 0.2e1 * t26, t24 ^ 2 + t26 ^ 2, t43 ^ 2, t43 * t52, 0, 0, 0, t23 * t52, -0.2e1 * t23 * t43, -0.2e1 * t3 * t21 - 0.2e1 * t47 * t5, t14 ^ 2 + t3 ^ 2 + t5 ^ 2, t37 * t18, -0.2e1 * t18 * t54, -0.2e1 * t47 * t53, 0.2e1 * t47 * t55, t17, 0.2e1 * t1 * t47 - 0.2e1 * t3 * t55, -0.2e1 * t2 * t47 - 0.2e1 * t3 * t53; 0, 0, 0, -1, 0, -pkin(1), -t41, t39, -t24 * t41 + t26 * t39, 0, 0, 0, 0, 0, -t41 * t56, t43 * t41, -t11 * t21 - t13 * t47, t3 * t11 + t5 * t13 - t14 * t41, 0, 0, 0, 0, 0, -t11 * t55 + t47 * t7, -t11 * t53 - t47 * t8; 0, 0, 0, 0, 0, 1, 0, 0, t39 ^ 2 + t35, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t13 ^ 2 + t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t21 - t3 * t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t47 + t13 * t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t17 + t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t56, 0, -t43 * t51, -t48 (t21 * t40 - t38 * t47) * pkin(4) (-t3 * t40 + t38 * t5) * pkin(4), -t42 * t53 (t36 - t37) * t21, t16, t15, 0, -t3 * t44 + t50 * t42, t3 * t42 + t50 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t39, -t56 * t39, 0 (-t11 * t40 + t13 * t38) * pkin(4), 0, 0, 0, 0, 0, -t11 * t44, t11 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t43, 0 (t21 * t38 + t40 * t47) * pkin(4), 0, 0, 0, 0, 0, t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t38 ^ 2 + t40 ^ 2) * pkin(4) ^ 2, t36, 0.2e1 * t54, 0, 0, 0, -0.2e1 * t29 * t44, 0.2e1 * t29 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t55, t47, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t44, 0, -t42 * t28, -t44 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t4;
