% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:53:44
% EndTime: 2019-05-05 18:53:46
% DurationCPUTime: 0.55s
% Computational Cost: add. (538->81), mult. (1015->131), div. (0->0), fcn. (1262->8), ass. (0->60)
t36 = sin(pkin(10));
t37 = cos(pkin(10));
t40 = sin(qJ(3));
t60 = cos(qJ(3));
t20 = t40 * t36 - t60 * t37;
t21 = t60 * t36 + t40 * t37;
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t12 = -t42 * t20 + t39 * t21;
t68 = 0.2e1 * t12;
t67 = -0.2e1 * t21;
t30 = -t37 * pkin(2) - pkin(1);
t66 = 0.2e1 * t30;
t38 = sin(qJ(6));
t65 = -0.2e1 * t38;
t41 = cos(qJ(6));
t64 = 0.2e1 * t41;
t52 = pkin(7) + qJ(2);
t26 = t52 * t36;
t27 = t52 * t37;
t14 = t60 * t26 + t40 * t27;
t8 = -t21 * pkin(8) + t14;
t15 = -t40 * t26 + t60 * t27;
t9 = t20 * pkin(8) + t15;
t4 = t39 * t9 - t42 * t8;
t63 = t4 * t38;
t62 = t4 * t41;
t43 = -pkin(3) - pkin(4);
t24 = t39 * qJ(4) - t42 * t43;
t22 = pkin(5) + t24;
t61 = pkin(5) + t22;
t59 = t38 * t12;
t13 = t39 * t20 + t42 * t21;
t58 = t38 * t13;
t57 = t38 * t41;
t56 = t41 * t12;
t55 = t41 * t13;
t54 = t42 * t38;
t53 = t42 * t41;
t51 = t36 ^ 2 + t37 ^ 2;
t50 = -0.2e1 * t13 * t12;
t49 = -0.2e1 * t57;
t48 = t38 * t55;
t10 = t20 * pkin(3) - t21 * qJ(4) + t30;
t47 = -pkin(5) * t13 - pkin(9) * t12;
t25 = t42 * qJ(4) + t39 * t43;
t23 = -pkin(9) + t25;
t46 = -t12 * t23 + t13 * t22;
t45 = -t12 * t39 - t13 * t42;
t7 = -t20 * pkin(4) - t10;
t35 = t41 ^ 2;
t34 = t38 ^ 2;
t28 = 0.2e1 * t57;
t11 = t13 ^ 2;
t6 = (t34 - t35) * t13;
t5 = t39 * t8 + t42 * t9;
t3 = t12 * pkin(5) - t13 * pkin(9) + t7;
t2 = t38 * t3 + t41 * t5;
t1 = t41 * t3 - t38 * t5;
t16 = [1, 0, 0, 0.2e1 * pkin(1) * t37, -0.2e1 * pkin(1) * t36, 0.2e1 * t51 * qJ(2), t51 * qJ(2) ^ 2 + pkin(1) ^ 2, t21 ^ 2, t20 * t67, 0, 0, 0, t20 * t66, t21 * t66, 0.2e1 * t10 * t20, 0.2e1 * t14 * t21 - 0.2e1 * t15 * t20, t10 * t67, t10 ^ 2 + t14 ^ 2 + t15 ^ 2, t11, t50, 0, 0, 0, t7 * t68, 0.2e1 * t7 * t13, t35 * t11, t11 * t49, t55 * t68, t38 * t50, t12 ^ 2, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t58, -0.2e1 * t2 * t12 + 0.2e1 * t4 * t55; 0, 0, 0, -t37, t36, 0, -pkin(1), 0, 0, 0, 0, 0, t20, t21, t20, 0, -t21, t10, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t56, t59; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t14, -t15, -t14, -pkin(3) * t21 - t20 * qJ(4), t15, -t14 * pkin(3) + t15 * qJ(4), 0, 0, -t13, t12, 0, t4, t5, -t48, t6, -t59, -t56, 0, t46 * t38 + t62, t46 * t41 - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t24, 0.2e1 * t25, t34, t28, 0, 0, 0, t22 * t64, t22 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t38, t45 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t42, t39, 0, 0, 0, 0, 0, -t53, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t4, -t5, t48, -t6, t59, t56, 0, t47 * t38 - t62, t47 * t41 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t24, -t25, -t34, t49, 0, 0, 0, -t61 * t41, t61 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t39, 0, 0, 0, 0, 0, t53, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t34, t28, 0, 0, 0, pkin(5) * t64, pkin(5) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t58, t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t41, 0, -t38 * t23, -t41 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t39, -t41 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * pkin(9), -t41 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t16;
