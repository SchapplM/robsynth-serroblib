% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:57
% EndTime: 2019-12-31 19:54:58
% DurationCPUTime: 0.49s
% Computational Cost: add. (1000->93), mult. (2242->163), div. (0->0), fcn. (2054->6), ass. (0->59)
t41 = sin(qJ(2));
t70 = sin(pkin(8));
t60 = t70 * t41;
t42 = cos(qJ(2));
t71 = cos(pkin(8));
t61 = t71 * t42;
t53 = t60 - t61;
t74 = t71 * t41 + t70 * t42;
t43 = 2 * qJD(5);
t73 = cos(qJ(4));
t72 = -qJ(3) - pkin(6);
t19 = t53 * t72;
t40 = sin(qJ(4));
t69 = qJD(4) * t40;
t68 = t41 * qJD(2);
t67 = t42 * qJD(2);
t66 = -0.2e1 * pkin(1) * qJD(2);
t39 = pkin(2) * t68;
t65 = -t42 * pkin(2) - pkin(1);
t64 = t70 * pkin(2);
t63 = qJD(4) * t73;
t58 = qJD(2) * t72;
t57 = t40 * t64;
t37 = qJD(2) * t60;
t56 = qJD(2) * t61 - t37;
t38 = t71 * pkin(2) + pkin(3);
t24 = qJD(4) * t57 - t38 * t63;
t55 = -t41 * qJD(3) + t42 * t58;
t54 = t42 * qJD(3) + t41 * t58;
t51 = t40 * t38 + t73 * t64;
t50 = t74 * qJD(2);
t49 = t73 * t53;
t20 = pkin(3) * t50 + t39;
t18 = t74 * t72;
t23 = t53 * pkin(3) + t65;
t17 = -t40 * t53 + t73 * t74;
t48 = -pkin(7) * t74 + t18;
t47 = t40 * t48;
t46 = t73 * t48;
t14 = t71 * t54 + t70 * t55;
t13 = -t70 * t54 + t71 * t55;
t45 = t56 * pkin(7) - t13;
t44 = -pkin(7) * t50 + t14;
t30 = -t73 * t38 - pkin(4) + t57;
t29 = qJ(5) + t51;
t25 = t51 * qJD(4);
t22 = qJD(5) - t24;
t21 = 0.2e1 * t25;
t16 = t40 * t74 + t49;
t15 = -t53 * pkin(7) + t19;
t8 = t17 * qJD(4) + t40 * t56 + t73 * t50;
t7 = qJD(4) * t49 + t40 * t50 - t73 * t56 + t69 * t74;
t6 = t73 * t15 + t47;
t5 = t40 * t15 - t46;
t4 = t16 * pkin(4) - t17 * qJ(5) + t23;
t3 = t8 * pkin(4) + t7 * qJ(5) - t17 * qJD(5) + t20;
t2 = qJD(4) * t47 + t15 * t63 + t40 * t44 + t73 * t45;
t1 = -qJD(4) * t46 + t15 * t69 + t40 * t45 - t73 * t44;
t9 = [0, 0, 0, 0.2e1 * t41 * t67, 0.2e1 * (-t41 ^ 2 + t42 ^ 2) * qJD(2), 0, 0, 0, t41 * t66, t42 * t66, -0.2e1 * t13 * t74 - 0.2e1 * t14 * t53 - 0.2e1 * t18 * t56 - 0.2e1 * t19 * t50, 0.2e1 * t18 * t13 + 0.2e1 * t19 * t14 + 0.2e1 * t65 * t39, -0.2e1 * t17 * t7, 0.2e1 * t7 * t16 - 0.2e1 * t17 * t8, 0, 0, 0, 0.2e1 * t20 * t16 + 0.2e1 * t23 * t8, 0.2e1 * t20 * t17 - 0.2e1 * t23 * t7, 0.2e1 * t3 * t16 + 0.2e1 * t4 * t8, 0.2e1 * t1 * t16 + 0.2e1 * t2 * t17 - 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t3 * t17 + 0.2e1 * t4 * t7, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, t67, -t68, 0, -pkin(6) * t67, pkin(6) * t68, (t71 * t37 + (-t71 ^ 2 * t42 - t70 * t74) * qJD(2)) * pkin(2), (t71 * t13 + t70 * t14) * pkin(2), 0, 0, -t7, -t8, 0, -t2, t1, -t2, -t22 * t16 + t25 * t17 - t29 * t8 - t30 * t7, -t1, -t1 * t29 + t2 * t30 + t6 * t22 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0.2e1 * t24, -t21, 0, 0.2e1 * t22, 0.2e1 * t29 * t22 + 0.2e1 * t30 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, t8, -t7, t8, 0, t7, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t2, t1, -t2, pkin(4) * t7 - t8 * qJ(5) - t16 * qJD(5), -t1, -t2 * pkin(4) - t1 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, -t25, 0, t43 - t24, -t25 * pkin(4) + t22 * qJ(5) + t29 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, qJ(5) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
