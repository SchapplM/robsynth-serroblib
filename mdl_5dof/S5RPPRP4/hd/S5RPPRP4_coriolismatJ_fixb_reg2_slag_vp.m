% Calculate inertial parameters regressor of coriolis matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:24
% EndTime: 2019-12-31 17:52:26
% DurationCPUTime: 0.81s
% Computational Cost: add. (490->88), mult. (767->100), div. (0->0), fcn. (634->4), ass. (0->73)
t57 = sin(qJ(4));
t53 = t57 ^ 2;
t58 = cos(qJ(4));
t54 = t58 ^ 2;
t41 = t53 + t54;
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t59 = -pkin(1) - pkin(2);
t91 = t56 * qJ(2) + t55 * t59;
t33 = -pkin(6) + t91;
t95 = -qJ(5) + t33;
t93 = pkin(4) * t57;
t92 = t58 * pkin(4);
t62 = -t55 * qJ(2) + t56 * t59;
t32 = pkin(3) - t62;
t23 = t32 + t92;
t1 = t23 * t93;
t90 = t1 * qJD(1);
t19 = t95 * t57;
t20 = t95 * t58;
t9 = t19 * t57 + t20 * t58;
t2 = t23 * t55 + t9 * t56;
t89 = t2 * qJD(1);
t3 = t56 * t93;
t88 = t3 * qJD(1);
t31 = t41 * t56;
t5 = t33 * t31 + t32 * t55;
t87 = t5 * qJD(1);
t86 = t9 * qJD(1);
t10 = -t62 * t55 + t91 * t56;
t85 = t10 * qJD(1);
t14 = (t23 + t92) * t57;
t84 = t14 * qJD(1);
t17 = t53 * pkin(4) - t23 * t58;
t83 = t17 * qJD(1);
t82 = t20 * qJD(4);
t21 = (0.1e1 / 0.2e1 + t54 / 0.2e1 + t53 / 0.2e1) * t55;
t81 = t21 * qJD(1);
t80 = t31 * qJD(1);
t79 = t41 * qJD(1);
t42 = t54 - t53;
t78 = t42 * qJD(1);
t77 = t55 * qJD(1);
t76 = t55 * qJD(2);
t75 = t56 * qJD(1);
t74 = t56 * qJD(2);
t73 = t57 * qJD(1);
t49 = t57 * qJD(4);
t72 = t58 * qJD(1);
t71 = t58 * qJD(4);
t70 = qJ(2) * qJD(1);
t69 = pkin(4) * t49;
t68 = pkin(4) * t73;
t67 = t57 * t76;
t66 = t56 * t72;
t65 = t55 * t73;
t64 = t56 * t73;
t63 = t55 * t71;
t61 = qJD(5) - t74;
t60 = qJD(1) * t32 - t74;
t44 = t57 * t71;
t43 = t57 * t72;
t40 = t55 * t72;
t39 = t58 * t76;
t34 = t42 * qJD(4);
t29 = t31 * qJD(2);
t27 = t56 * t49 - t40;
t26 = t55 * t49 + t66;
t25 = -t63 + t64;
t24 = t56 * t71 + t65;
t22 = (0.1e1 / 0.2e1 - t41 / 0.2e1) * t55;
t13 = (-0.1e1 + t41) * t55 * t74;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, 0, 0, 0, 0, 0, t76, t74, 0, t10 * qJD(2), t44, t34, 0, -t44, 0, 0, -t32 * t49 + t39, -t32 * t71 - t67, -t29, t5 * qJD(2), t44, t34, 0, -t44, 0, 0, -t14 * qJD(4) + t39, t17 * qJD(4) - t67, t41 * qJD(5) - t29, t2 * qJD(2) - t1 * qJD(4) - t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t70, 0, 0, 0, 0, 0, 0, t77, t75, 0, t85, 0, 0, 0, 0, 0, 0, t40, -t65, -t80, t13 + t87, 0, 0, 0, 0, 0, 0, t40, -t65, -t80, t22 * qJD(5) + t13 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t78, -t71, -t43, t49, 0, -t32 * t73 - t33 * t71, -t32 * t72 + t33 * t49, 0, 0, t43, t78, -t71, -t43, t49, 0, -t82 - t84, t19 * qJD(4) + t83, pkin(4) * t71, -pkin(4) * t82 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t22 * qJD(2) - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t70, 0, 0, 0, 0, 0, 0, -t77, -t75, 0, -t85, 0, 0, 0, 0, 0, 0, t27, t24, t80, -t87, 0, 0, 0, 0, 0, 0, t27, t24, t80, t3 * qJD(4) - t21 * qJD(5) - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, 0, -pkin(4) * t63 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t71, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t71, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t78, 0, t43, 0, 0, t60 * t57, t60 * t58, 0, 0, -t43, -t78, 0, t43, 0, 0, t61 * t57 + t84, t61 * t58 - t83, 0, -t3 * qJD(2) + qJD(5) * t93 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t66, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t66, 0, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t71, -t79, t21 * qJD(2) - t69 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t72, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
