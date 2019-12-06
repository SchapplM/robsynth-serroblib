% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:59
% EndTime: 2019-12-05 18:41:01
% DurationCPUTime: 0.62s
% Computational Cost: add. (657->79), mult. (1826->135), div. (0->0), fcn. (1258->8), ass. (0->67)
t56 = cos(qJ(3));
t57 = cos(qJ(2));
t68 = t57 * pkin(1) + pkin(2);
t38 = t56 * t68;
t77 = pkin(1) * qJD(2);
t72 = t57 * t77;
t53 = sin(qJ(3));
t54 = sin(qJ(2));
t74 = t53 * t54 * pkin(1);
t84 = qJD(2) + qJD(3);
t15 = -qJD(3) * t38 - t56 * t72 + t84 * t74;
t51 = cos(pkin(9));
t80 = t56 * t54;
t58 = t84 * (-t53 * t57 - t80) * pkin(1);
t76 = qJD(3) * pkin(2);
t71 = t53 * t76;
t16 = t58 - t71;
t50 = sin(pkin(9));
t82 = t50 * t16;
t6 = -t51 * t15 + t82;
t52 = sin(qJ(5));
t48 = t52 ^ 2;
t55 = cos(qJ(5));
t49 = t55 ^ 2;
t78 = t48 + t49;
t87 = t78 * t6;
t39 = t50 * t71;
t70 = t56 * t76;
t27 = t51 * t70 - t39;
t86 = t78 * t27;
t81 = t51 * t53;
t26 = (t50 * t56 + t81) * t76;
t5 = -t50 * t15 - t51 * t16;
t85 = -t26 - t5;
t47 = t55 * qJD(5);
t31 = t38 - t74;
t28 = pkin(3) + t31;
t32 = pkin(1) * t80 + t53 * t68;
t12 = t51 * t28 - t50 * t32;
t9 = -pkin(4) - t12;
t83 = t9 * t47 + t5 * t52;
t46 = t56 * pkin(2) + pkin(3);
t29 = -t50 * t53 * pkin(2) + t51 * t46;
t24 = -pkin(4) - t29;
t79 = t24 * t47 + t26 * t52;
t13 = t50 * t28 + t51 * t32;
t30 = pkin(2) * t81 + t50 * t46;
t75 = t52 * qJD(5);
t73 = t54 * t77;
t69 = t52 * t47;
t7 = t9 * t75;
t67 = -t5 * t55 + t7;
t44 = t50 * pkin(3) + pkin(8);
t64 = t78 * t44;
t17 = t24 * t75;
t63 = -t26 * t55 + t17;
t62 = -0.2e1 * t71;
t61 = t15 - t70;
t45 = -t51 * pkin(3) - pkin(4);
t42 = -0.2e1 * t69;
t41 = 0.2e1 * t69;
t35 = t45 * t47;
t34 = t45 * t75;
t33 = 0.2e1 * (-t48 + t49) * qJD(5);
t25 = pkin(8) + t30;
t10 = pkin(8) + t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t73, -0.2e1 * t72, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16, 0.2e1 * t15, 0, -0.2e1 * t32 * t15 + 0.2e1 * t31 * t16, 0, 0, 0, 0, 0, 0, -0.2e1 * t5, -0.2e1 * t6, 0, -0.2e1 * t12 * t5 + 0.2e1 * t13 * t6, t41, t33, 0, t42, 0, 0, 0.2e1 * t67, 0.2e1 * t83, 0.2e1 * t87, 0.2e1 * t10 * t87 + 0.2e1 * t9 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t72, 0, 0, 0, 0, 0, 0, 0, 0, t58 + t62, t61, 0, (-t15 * t53 + t16 * t56 + (-t31 * t53 + t32 * t56) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t85, t61 * t51 + t39 - t82, 0, -t12 * t26 + t13 * t27 - t5 * t29 + t6 * t30, t41, t33, 0, t42, 0, 0, t85 * t55 + t17 + t7, t79 + t83, t86 + t87, t10 * t86 + t5 * t24 + t25 * t87 + t9 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -0.2e1 * t70, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26, -0.2e1 * t27, 0, -0.2e1 * t29 * t26 + 0.2e1 * t30 * t27, t41, t33, 0, t42, 0, 0, 0.2e1 * t63, 0.2e1 * t79, 0.2e1 * t86, 0.2e1 * t24 * t26 + 0.2e1 * t25 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, (-t5 * t51 + t50 * t6) * pkin(3), t41, t33, 0, t42, 0, 0, t34 + t67, t35 + t83, t87, t5 * t45 + t6 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, (-t26 * t51 + t27 * t50) * pkin(3), t41, t33, 0, t42, 0, 0, t34 + t63, t35 + t79, t86, t26 * t45 + t27 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t33, 0, t42, 0, 0, 0.2e1 * t34, 0.2e1 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t75, 0, -t10 * t47 - t52 * t6, t10 * t75 - t55 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t75, 0, -t25 * t47 - t52 * t27, t25 * t75 - t55 * t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t75, 0, -t44 * t47, t44 * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
