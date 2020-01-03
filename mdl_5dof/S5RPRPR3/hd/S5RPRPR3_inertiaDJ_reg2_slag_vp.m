% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:22
% EndTime: 2020-01-03 11:36:26
% DurationCPUTime: 0.70s
% Computational Cost: add. (690->90), mult. (1495->164), div. (0->0), fcn. (1132->8), ass. (0->76)
t44 = cos(pkin(8)) * pkin(1) + pkin(2);
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t97 = sin(pkin(8)) * pkin(1);
t100 = -t55 * t44 + t53 * t97;
t80 = t53 * t44 + t55 * t97;
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t99 = t51 * pkin(4) + t49 * pkin(7);
t98 = 2 * qJD(5);
t25 = t100 * qJD(3);
t24 = qJD(4) - t25;
t26 = t80 * qJD(3);
t28 = qJ(4) + t80;
t52 = sin(qJ(5));
t74 = qJD(5) * t52;
t40 = t51 * t74;
t54 = cos(qJ(5));
t59 = -pkin(3) + t100;
t57 = t59 - t99;
t56 = t54 * t57;
t87 = t51 * t54;
t2 = -qJD(5) * t56 - t24 * t87 - t52 * t26 + t28 * t40;
t96 = t2 * t52;
t61 = pkin(3) + t99;
t58 = t54 * t61;
t64 = qJ(4) * t74;
t72 = t51 * qJD(4);
t9 = qJD(5) * t58 + t51 * t64 - t54 * t72;
t93 = t9 * t52;
t71 = t52 * qJD(4);
t77 = qJ(4) * t51;
t23 = -t52 * t61 + t54 * t77;
t75 = qJD(5) * t23;
t10 = -t51 * t71 - t75;
t6 = t28 * t87 + t52 * t57;
t76 = qJD(5) * t6;
t86 = t52 * t24;
t3 = t54 * t26 - t51 * t86 - t76;
t92 = -t10 - t3;
t47 = t49 ^ 2;
t18 = t47 * t24;
t91 = t54 * t18 - t2 * t51;
t45 = t47 * qJD(4);
t90 = t54 * t45 - t9 * t51;
t89 = t26 * t49;
t88 = t26 * t51;
t48 = t51 ^ 2;
t19 = t48 * t24;
t73 = qJD(5) * t54;
t67 = t47 * t73;
t84 = t28 * t67 + t47 * t86;
t78 = qJ(4) * t24;
t83 = t28 * t45 + t47 * t78;
t82 = t19 + t18;
t81 = qJ(4) * t67 + t47 * t71;
t79 = t48 * qJD(4) + t45;
t70 = qJ(4) * qJD(4);
t68 = t47 * t74;
t66 = t49 * t74;
t65 = t49 * t73;
t63 = t49 * t51 * t98;
t62 = t52 * t67;
t43 = t47 * t70;
t41 = t51 * t73;
t33 = -0.2e1 * t62;
t32 = 0.2e1 * t62;
t31 = t54 * t63;
t30 = t52 * t63;
t27 = (t52 ^ 2 - t54 ^ 2) * t47 * t98;
t22 = -t52 * t77 - t58;
t12 = t22 * t66;
t7 = t28 * t18;
t5 = -t52 * t51 * t28 + t56;
t4 = t5 * t66;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26, 0.2e1 * t25, 0, 0.2e1 * t100 * t26 - 0.2e1 * t80 * t25, 0, 0, 0, 0, 0, 0, -0.2e1 * t88, 0.2e1 * t89, 0.2e1 * t82, 0.2e1 * t28 * t19 + 0.2e1 * t59 * t26 + 0.2e1 * t7, t33, t27, t30, t32, t31, 0, -0.2e1 * t3 * t51 + 0.2e1 * t84, -0.2e1 * t28 * t68 + 0.2e1 * t91, 0.2e1 * t4 + 0.2e1 * (t96 + (-t3 - t76) * t54) * t49, -0.2e1 * t6 * t2 + 0.2e1 * t5 * t3 + 0.2e1 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t2 * t54 - t51 * t24 - t3 * t52 + (-t5 * t54 - t52 * t6) * qJD(5)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t89, t79 + t82, -t26 * pkin(3) + (qJD(4) * t28 + t78) * t48 + t83, t33, t27, t30, t32, t31, 0, t92 * t51 + t81 + t84, (-qJ(4) - t28) * t68 + t90 + t91, t12 + t4 + ((t2 + t9) * t52 + ((-t23 - t6) * qJD(5) + t92) * t54) * t49, t5 * t10 - t2 * t23 + t3 * t22 - t6 * t9 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t72 - t10 * t52 - t54 * t9 + (-t22 * t54 - t23 * t52) * qJD(5)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79, 0.2e1 * t48 * t70 + 0.2e1 * t43, t33, t27, t30, t32, t31, 0, -0.2e1 * t10 * t51 + 0.2e1 * t81, -0.2e1 * t47 * t64 + 0.2e1 * t90, 0.2e1 * t12 + 0.2e1 * (t93 + (-t10 - t75) * t54) * t49, 0.2e1 * t22 * t10 - 0.2e1 * t23 * t9 + 0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, t40, t41, 0, -t96 + t3 * t54 + (-t5 * t52 + t54 * t6) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t41, 0, t10 * t54 - t93 + (-t22 * t52 + t23 * t54) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, -t65, 0, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t66, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, -t65, 0, t10, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
