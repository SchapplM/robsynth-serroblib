% Calculate inertial parameters regressor of coriolis matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:43
% EndTime: 2019-12-31 16:41:44
% DurationCPUTime: 0.38s
% Computational Cost: add. (594->54), mult. (1010->66), div. (0->0), fcn. (1083->4), ass. (0->48)
t50 = sin(pkin(6));
t46 = t50 ^ 2;
t51 = cos(pkin(6));
t47 = t51 ^ 2;
t41 = t46 + t47;
t72 = sin(qJ(4));
t73 = cos(qJ(4));
t34 = t73 * t50 + t72 * t51;
t32 = t34 ^ 2;
t36 = -t72 * t50 + t73 * t51;
t74 = t36 ^ 2;
t71 = -pkin(1) - qJ(3);
t58 = -pkin(5) + t71;
t38 = t58 * t50;
t53 = t58 * t51;
t16 = t72 * t38 - t73 * t53;
t17 = t73 * t38 + t72 * t53;
t7 = t16 * t36 - t17 * t34;
t70 = t7 * qJD(1);
t52 = -t32 / 0.2e1 - t74 / 0.2e1;
t9 = -0.1e1 / 0.2e1 + t52;
t69 = t9 * qJD(1);
t15 = t32 - t74;
t68 = t15 * qJD(1);
t18 = t32 + t74;
t67 = t18 * qJD(1);
t33 = t41 * t71;
t66 = t33 * qJD(1);
t65 = t34 * qJD(1);
t26 = t34 * qJD(4);
t64 = t36 * qJD(1);
t29 = t36 * qJD(4);
t55 = -t46 / 0.2e1 - t47 / 0.2e1;
t40 = -0.1e1 / 0.2e1 + t55;
t63 = t40 * qJD(1);
t62 = t41 * qJD(1);
t42 = t50 * pkin(3) + qJ(2);
t61 = t42 * qJD(1);
t60 = t50 * qJD(1);
t59 = t51 * qJD(1);
t57 = t34 * t64;
t56 = t34 * t29;
t54 = qJD(3) + t61;
t49 = qJ(2) * qJD(2);
t48 = qJD(1) * qJ(2);
t39 = 0.1e1 / 0.2e1 + t55;
t8 = 0.1e1 / 0.2e1 + t52;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t49, 0, 0, 0, 0, 0, 0, qJD(2) * t50, qJD(2) * t51, t41 * qJD(3), -t33 * qJD(3) + t49, -t56, t15 * qJD(4), 0, t56, 0, 0, qJD(2) * t34 + t42 * t29, qJD(2) * t36 - t42 * t26, t18 * qJD(3), t42 * qJD(2) + t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t48, 0, 0, 0, 0, 0, 0, t60, t59, 0, t39 * qJD(3) + t48, 0, 0, 0, 0, 0, 0, t65, t64, 0, t8 * qJD(3) + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t39 * qJD(2) - t66, 0, 0, 0, 0, 0, 0, 0, 0, t67, t8 * qJD(2) + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t68, -t26, t57, -t29, 0, -t17 * qJD(4) + t36 * t61, t16 * qJD(4) - t34 * t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t48, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, t40 * qJD(3) - t48, 0, 0, 0, 0, 0, 0, -t65, -t64, 0, t9 * qJD(3) - t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t40 * qJD(2) + t66, 0, 0, 0, 0, 0, 0, t29, -t26, -t67, -t9 * qJD(2) - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t68, 0, -t57, 0, 0, -t54 * t36, t54 * t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
