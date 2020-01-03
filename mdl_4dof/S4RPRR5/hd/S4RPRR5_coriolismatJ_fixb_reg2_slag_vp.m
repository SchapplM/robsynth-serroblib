% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:42
% EndTime: 2019-12-31 16:51:43
% DurationCPUTime: 0.59s
% Computational Cost: add. (588->100), mult. (1001->134), div. (0->0), fcn. (724->4), ass. (0->72)
t54 = sin(qJ(4));
t52 = t54 ^ 2;
t56 = cos(qJ(4));
t53 = t56 ^ 2;
t91 = t52 + t53;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t97 = -pkin(1) - pkin(2);
t38 = t57 * qJ(2) + t55 * t97;
t72 = qJD(1) - qJD(3);
t105 = t72 * t38;
t45 = t53 - t52;
t104 = t72 * t45;
t36 = -pkin(6) + t38;
t37 = t55 * qJ(2) - t57 * t97;
t66 = t53 / 0.2e1 + t52 / 0.2e1;
t64 = t66 * pkin(6);
t35 = pkin(3) + t37;
t71 = -pkin(3) / 0.2e1 - t35 / 0.2e1;
t1 = (t66 * t37 + t71) * t55 + (t38 / 0.2e1 - t66 * t36 + t64) * t57;
t15 = (-0.1e1 + t91) * t57 * t55;
t14 = t15 * qJD(2);
t103 = t1 * qJD(1) - t14;
t102 = t35 + t37;
t34 = t91 * t57;
t83 = t34 * qJD(1);
t100 = -t34 * qJD(3) + t83;
t96 = pkin(3) * t54;
t95 = t55 * pkin(3);
t94 = t35 * t55;
t93 = t37 * t55;
t92 = t38 * t57;
t12 = t91 * t37;
t3 = -t36 * t12 + t35 * t38;
t89 = t3 * qJD(1);
t4 = t36 * t34 + t94;
t88 = t4 * qJD(1);
t87 = qJD(1) * t35;
t86 = qJD(1) * t56;
t85 = qJD(3) * t56;
t13 = t92 + t93;
t84 = t13 * qJD(1);
t82 = t34 * qJD(2);
t80 = t45 * qJD(4);
t79 = t54 * qJD(4);
t78 = t55 * qJD(1);
t77 = t55 * qJD(2);
t76 = t56 * qJD(4);
t75 = t57 * qJD(1);
t74 = t57 * qJD(2);
t73 = qJ(2) * qJD(1);
t70 = t54 * t78;
t69 = t56 * t78;
t39 = t72 * t55;
t65 = t72 * t56;
t40 = t72 * t57;
t63 = t37 / 0.2e1 + t71;
t62 = -t74 + t87;
t61 = t38 * qJD(1) + t77;
t60 = t38 * qJD(3) + t77;
t6 = t63 * t54;
t59 = t6 * qJD(1) + qJD(3) * t96;
t7 = t63 * t56;
t58 = pkin(3) * t85 + t7 * qJD(1);
t47 = t54 * t76;
t31 = (-t85 + t86) * t54;
t17 = t55 * t65 - t57 * t79;
t16 = t54 * t39 + t57 * t76;
t9 = (pkin(3) + t102) * t56 / 0.2e1;
t8 = t96 / 0.2e1 + t102 * t54 / 0.2e1;
t2 = -t92 / 0.2e1 + t94 / 0.2e1 - t95 / 0.2e1 + t57 * t64 + t91 * (-t93 / 0.2e1 + t36 * t57 / 0.2e1);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, 0, 0, 0, 0, 0, t60, -t37 * qJD(3) + t74, 0, t13 * qJD(2), t47, t80, 0, -t47, 0, 0, -t35 * t79 + t60 * t56, -t35 * t76 - t60 * t54, t12 * qJD(3) - t82, t4 * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t73, 0, 0, 0, 0, 0, 0, t78, t75, 0, t84, 0, 0, 0, 0, 0, 0, t69, -t70, -t83, t2 * qJD(3) + t14 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t72 * t37, 0, 0, -t47, -t80, 0, t47, 0, 0, t8 * qJD(4) + t38 * t65, t9 * qJD(4) - t54 * t105, t72 * t12, t89 + t2 * qJD(2) + (-t38 * pkin(3) - pkin(6) * t12) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t104, -t76, -t31, t79, 0, t8 * qJD(3) - t36 * t76 - t54 * t87, t9 * qJD(3) - t35 * t86 + t36 * t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t73, 0, 0, 0, 0, 0, 0, -t39, -t40, 0, -t84, 0, 0, 0, 0, 0, 0, -t17, t16, t100, -t1 * qJD(3) - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t40, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t100, (pkin(6) * t34 - t95) * qJD(3) - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t40 - t55 * t76, t56 * t40 + t55 * t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t37 * qJD(1) - t74, 0, 0, -t47, -t80, 0, t47, 0, 0, -t6 * qJD(4) - t61 * t56, -t7 * qJD(4) + t61 * t54, -t12 * qJD(1) + t82, t1 * qJD(2) - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t75, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t70, t83, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t80, 0, -t47, 0, 0, -pkin(3) * t79, -pkin(3) * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t104, t76, t31, -t79, 0, -pkin(6) * t76 - t59, pkin(6) * t79 - t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t104, 0, t31, 0, 0, t6 * qJD(3) + t62 * t54, t7 * qJD(3) + t62 * t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * t75, -t56 * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t104, 0, -t31, 0, 0, t59, t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
