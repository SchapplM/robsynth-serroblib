% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:44
% DurationCPUTime: 0.83s
% Computational Cost: add. (373->104), mult. (962->160), div. (0->0), fcn. (563->4), ass. (0->71)
t60 = sin(qJ(3));
t57 = t60 ^ 2;
t62 = cos(qJ(3));
t58 = t62 ^ 2;
t63 = cos(qJ(2));
t84 = pkin(1) * qJD(2);
t76 = t63 * t84;
t96 = (t57 + t58) * t76;
t53 = t62 * qJD(4);
t83 = t60 * qJ(4);
t68 = -t62 * pkin(3) - t83;
t94 = t68 * qJD(3) + t53;
t59 = -pkin(3) - qJ(5);
t93 = t59 * t62 - t83;
t92 = 0.2e1 * (-t57 + t58) * qJD(3);
t52 = t60 * qJD(3);
t72 = pkin(3) * t52 - t60 * qJD(4);
t81 = qJD(3) * qJ(4);
t4 = qJ(5) * t52 + (-qJD(5) - t81) * t62 + t72;
t61 = sin(qJ(2));
t51 = t61 * t84;
t1 = t51 + t4;
t91 = -t1 - t4;
t90 = t63 * pkin(1);
t47 = t61 * pkin(1) + pkin(7);
t89 = pkin(4) + t47;
t87 = t96 * t47;
t86 = t96 * pkin(7);
t48 = -pkin(2) - t90;
t54 = t62 * qJD(3);
t85 = t48 * t54 + t60 * t51;
t82 = qJ(4) * qJD(4);
t55 = t60 * pkin(4);
t25 = t60 * t47 + t55;
t42 = t62 * t76;
t6 = -t89 * t52 + t42;
t69 = t60 * t76;
t7 = t89 * t54 + t69;
t80 = t25 * t54 + t6 * t62 + t7 * t60;
t32 = (-pkin(4) - pkin(7)) * t52;
t49 = pkin(7) * t54;
t33 = pkin(4) * t54 + t49;
t38 = t60 * pkin(7) + t55;
t79 = t32 * t62 + t33 * t60 + t38 * t54;
t78 = pkin(2) * t52;
t77 = pkin(2) * t54;
t75 = pkin(7) * t52;
t74 = t60 * t54;
t36 = -pkin(2) + t68;
t24 = t36 - t90;
t73 = qJD(3) * (-t24 - t36);
t67 = t48 * t52 - t62 * t51;
t23 = -pkin(2) + t93;
t20 = -t62 * t81 + t72;
t64 = 0.2e1 * qJD(4);
t56 = t62 * pkin(4);
t45 = -0.2e1 * t74;
t44 = 0.2e1 * t74;
t39 = t62 * pkin(7) + t56;
t26 = t62 * t47 + t56;
t16 = t23 - t90;
t15 = t47 * t54 + t69;
t14 = t47 * t52 - t42;
t13 = t23 * t52;
t12 = 0.2e1 * t96;
t11 = t20 * t62;
t10 = t20 + t51;
t9 = qJD(3) * t93 - qJD(5) * t60 + t53;
t8 = t16 * t52;
t5 = t10 * t62;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t51, -0.2e1 * t76, 0, 0, t44, t92, 0, t45, 0, 0, 0.2e1 * t67, 0.2e1 * t85, t12, 0.2e1 * t48 * t51 + 0.2e1 * t87, 0, 0, 0, t44, t92, t45, t12, -0.2e1 * t24 * t52 + 0.2e1 * t5, -0.2e1 * t10 * t60 - 0.2e1 * t24 * t54, 0.2e1 * t24 * t10 + 0.2e1 * t87, 0, 0, 0, t45, -t92, t44, -0.2e1 * t26 * t52 + 0.2e1 * t80, -0.2e1 * t1 * t60 - 0.2e1 * t16 * t54, -0.2e1 * t1 * t62 + 0.2e1 * t8, 0.2e1 * t16 * t1 + 0.2e1 * t25 * t7 + 0.2e1 * t26 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t76, 0, 0, t44, t92, 0, t45, 0, 0, t67 - t78, -t77 + t85, t96, -pkin(2) * t51 + t86, 0, 0, 0, t44, t92, t45, t96, t60 * t73 + t11 + t5, (-t10 - t20) * t60 + t62 * t73, t10 * t36 + t24 * t20 + t86, 0, 0, 0, t45, -t92, t44, (-t26 - t39) * t52 + t79 + t80, t91 * t60 + (-t16 - t23) * t54, t91 * t62 + t13 + t8, t1 * t23 + t16 * t4 + t25 * t33 + t26 * t32 + t7 * t38 + t6 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t92, 0, t45, 0, 0, -0.2e1 * t78, -0.2e1 * t77, 0, 0, 0, 0, 0, t44, t92, t45, 0, -0.2e1 * t36 * t52 + 0.2e1 * t11, -0.2e1 * t20 * t60 - 0.2e1 * t36 * t54, 0.2e1 * t36 * t20, 0, 0, 0, t45, -t92, t44, -0.2e1 * t39 * t52 + 0.2e1 * t79, -0.2e1 * t23 * t54 - 0.2e1 * t4 * t60, -0.2e1 * t4 * t62 + 0.2e1 * t13, 0.2e1 * t23 * t4 + 0.2e1 * t39 * t32 + 0.2e1 * t38 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t52, 0, -t15, t14, 0, 0, 0, -t54, t52, 0, 0, 0, t94, t15, -t14, (-pkin(3) * t60 + qJ(4) * t62) * t76 + t94 * t47, 0, t52, t54, 0, 0, 0, t9, t6, -t7, t6 * qJ(4) + t26 * qJD(4) - t25 * qJD(5) + t7 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t52, 0, -t49, t75, 0, 0, 0, -t54, t52, 0, 0, 0, t94, t49, -t75, t94 * pkin(7), 0, t52, t54, 0, 0, 0, t9, t32, -t33, t32 * qJ(4) + t39 * qJD(4) - t38 * qJD(5) + t33 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0.2e1 * t82, 0, 0, 0, 0, 0, 0, 0, t64, 0.2e1 * qJD(5), -0.2e1 * t59 * qJD(5) + 0.2e1 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, t15, 0, 0, 0, 0, 0, 0, t54, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, t49, 0, 0, 0, 0, 0, 0, t54, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
