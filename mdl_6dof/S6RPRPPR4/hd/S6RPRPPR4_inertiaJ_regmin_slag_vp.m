% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:49:30
% EndTime: 2019-05-05 16:49:32
% DurationCPUTime: 0.61s
% Computational Cost: add. (700->91), mult. (1362->164), div. (0->0), fcn. (1591->8), ass. (0->66)
t52 = sin(pkin(9));
t54 = cos(pkin(9));
t56 = sin(qJ(3));
t77 = cos(qJ(3));
t30 = t56 * t52 - t77 * t54;
t81 = -0.2e1 * t30;
t51 = sin(pkin(10));
t53 = cos(pkin(10));
t40 = t51 ^ 2 + t53 ^ 2;
t33 = t77 * t52 + t56 * t54;
t83 = -0.2e1 * t33;
t55 = sin(qJ(6));
t57 = cos(qJ(6));
t32 = t57 * t51 - t55 * t53;
t66 = t51 * qJ(5) + pkin(3);
t78 = pkin(4) + pkin(5);
t26 = t78 * t53 + t66;
t82 = 0.2e1 * t26;
t43 = -t54 * pkin(2) - pkin(1);
t80 = 0.2e1 * t43;
t79 = -0.2e1 * t51;
t76 = t32 * t30;
t75 = t51 * t30;
t24 = t51 * t33;
t74 = t53 * t30;
t25 = t53 * t33;
t71 = pkin(7) + qJ(2);
t14 = t30 * pkin(3) - t33 * qJ(4) + t43;
t37 = t71 * t52;
t39 = t71 * t54;
t21 = -t56 * t37 + t77 * t39;
t9 = t51 * t14 + t53 * t21;
t70 = t40 * qJ(4) ^ 2;
t69 = t52 ^ 2 + t54 ^ 2;
t68 = qJ(4) * t30;
t67 = qJ(5) * t53;
t5 = t30 * qJ(5) + t9;
t16 = t51 * t21;
t8 = t53 * t14 - t16;
t19 = t77 * t37 + t56 * t39;
t6 = -t30 * pkin(4) - t8;
t65 = t5 * t53 + t6 * t51;
t64 = t5 * t51 - t6 * t53;
t63 = t9 * t51 + t8 * t53;
t62 = -t8 * t51 + t9 * t53;
t61 = -pkin(3) * t33 - t68;
t29 = t55 * t51 + t57 * t53;
t35 = -t53 * pkin(4) - t66;
t60 = -t33 * t35 + t68;
t45 = t51 * qJ(4);
t38 = (-pkin(8) + qJ(4)) * t53;
t36 = -t51 * pkin(8) + t45;
t34 = 0.2e1 * t40 * qJ(4);
t22 = t29 * t30;
t20 = t55 * t36 + t57 * t38;
t18 = t57 * t36 - t55 * t38;
t15 = t40 * t33;
t12 = t29 * t33;
t11 = t32 * t33;
t10 = (pkin(4) * t51 - t67) * t33 + t19;
t7 = (-t78 * t51 + t67) * t33 - t19;
t4 = pkin(8) * t24 + t5;
t3 = t16 + (-pkin(8) * t33 - t14) * t53 - t78 * t30;
t2 = t55 * t3 + t57 * t4;
t1 = t57 * t3 - t55 * t4;
t13 = [1, 0, 0, 0.2e1 * pkin(1) * t54, -0.2e1 * pkin(1) * t52, 0.2e1 * t69 * qJ(2), t69 * qJ(2) ^ 2 + pkin(1) ^ 2, t33 ^ 2, t33 * t81, 0, 0, 0, t30 * t80, t33 * t80, 0.2e1 * t19 * t24 + 0.2e1 * t8 * t30, 0.2e1 * t19 * t25 - 0.2e1 * t9 * t30, t63 * t83, t19 ^ 2 + t8 ^ 2 + t9 ^ 2, 0.2e1 * t10 * t24 - 0.2e1 * t6 * t30, t64 * t83, -0.2e1 * t10 * t25 + 0.2e1 * t5 * t30, t10 ^ 2 + t5 ^ 2 + t6 ^ 2, t12 ^ 2, 0.2e1 * t12 * t11, t12 * t81, t11 * t81, t30 ^ 2, -0.2e1 * t1 * t30 - 0.2e1 * t7 * t11, 0.2e1 * t7 * t12 + 0.2e1 * t2 * t30; 0, 0, 0, -t54, t52, 0, -pkin(1), 0, 0, 0, 0, 0, t30, t33, t74, -t75, -t15, t63, t74, -t15, t75, t64, 0, 0, 0, 0, 0, t22, t76; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30, 0, -t19, -t21, -t19 * t53 + t61 * t51, t19 * t51 + t61 * t53, t62, -t19 * pkin(3) + t62 * qJ(4), -t10 * t53 - t60 * t51, t65, -t10 * t51 + t60 * t53, t65 * qJ(4) + t10 * t35, t12 * t32, t32 * t11 - t12 * t29, -t76, t22, 0, -t26 * t11 - t18 * t30 + t7 * t29, t26 * t12 + t20 * t30 + t7 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t53, pkin(3) * t79, t34, pkin(3) ^ 2 + t70, -0.2e1 * t35 * t53, t34, t35 * t79, t35 ^ 2 + t70, t32 ^ 2, -0.2e1 * t32 * t29, 0, 0, 0, t29 * t82, t32 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, 0, t19, t24, 0, -t25, t10, 0, 0, 0, 0, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t51, 0, -pkin(3), -t53, 0, -t51, t35, 0, 0, 0, 0, 0, -t29, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t25, 0, t6, 0, 0, 0, 0, 0, -t57 * t30, t55 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, -t30, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, t18, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
