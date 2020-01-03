% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:27
% DurationCPUTime: 0.52s
% Computational Cost: add. (502->92), mult. (1261->133), div. (0->0), fcn. (804->6), ass. (0->67)
t42 = sin(qJ(4));
t40 = t42 ^ 2;
t45 = cos(qJ(4));
t41 = t45 ^ 2;
t68 = t40 + t41;
t43 = sin(qJ(3));
t44 = sin(qJ(2));
t77 = pkin(1) * t44;
t35 = t43 * t77;
t46 = cos(qJ(3));
t72 = cos(qJ(2));
t58 = t72 * pkin(1);
t53 = t58 + pkin(2);
t48 = t46 * t53;
t54 = qJD(2) * t58;
t79 = qJD(2) + qJD(3);
t8 = -qJD(3) * t48 + t79 * t35 - t46 * t54;
t57 = t68 * t8;
t67 = qJD(3) * pkin(2);
t61 = t46 * t67;
t82 = t68 * t61;
t71 = t46 * t44;
t81 = t79 * pkin(1) * (t43 * t72 + t71);
t78 = 2 * qJD(5);
t76 = t46 * pkin(2);
t39 = t45 * qJD(4);
t66 = t42 * qJD(4);
t22 = pkin(4) * t66 - qJ(5) * t39 - t42 * qJD(5);
t62 = t43 * t67;
t12 = t22 + t62;
t9 = t62 + t81;
t4 = t22 + t9;
t75 = -t12 - t4;
t74 = -t22 - t4;
t18 = t35 - t48 - pkin(3);
t73 = t18 * t39 + t9 * t42;
t70 = -t12 - t22;
t37 = -pkin(3) - t76;
t69 = t37 * t39 + t42 * t62;
t65 = pkin(3) * t66;
t64 = pkin(3) * t39;
t63 = qJD(2) * t77;
t60 = pkin(8) * t66;
t59 = pkin(8) * t39;
t13 = t18 * t66;
t56 = -t9 * t45 + t13;
t55 = -0.2e1 * t62;
t52 = -t45 * pkin(4) - t42 * qJ(5);
t51 = pkin(4) * t42 - qJ(5) * t45;
t26 = t37 * t66;
t49 = -t45 * t62 + t26;
t29 = -pkin(3) + t52;
t21 = t52 * qJD(4) + t45 * qJD(5);
t36 = t43 * pkin(2) + pkin(8);
t34 = 0.2e1 * t42 * t39;
t25 = 0.2e1 * (-t40 + t41) * qJD(4);
t24 = t29 - t76;
t23 = t29 * t66;
t19 = pkin(1) * t71 + t43 * t53 + pkin(8);
t17 = t24 * t66;
t16 = t36 * t39 + t42 * t61;
t15 = t36 * t66 - t45 * t61;
t11 = t18 + t52;
t10 = t11 * t66;
t3 = t19 * t39 - t42 * t8;
t2 = t19 * t66 + t45 * t8;
t1 = [0, 0, 0, 0, -0.2e1 * t63, -0.2e1 * t54, 0, -0.2e1 * t9, 0.2e1 * t8, t34, t25, 0, 0, 0, 0.2e1 * t56, 0.2e1 * t73, -0.2e1 * t4 * t45 + 0.2e1 * t10, -0.2e1 * t57, -0.2e1 * t11 * t39 - 0.2e1 * t4 * t42, 0.2e1 * t11 * t4 - 0.2e1 * t19 * t57; 0, 0, 0, 0, -t63, -t54, 0, t55 - t81, t8 - t61, t34, t25, 0, 0, 0, t13 + t26 + (-t9 - t62) * t45, t69 + t73, t75 * t45 + t10 + t17, t82 - t57, t75 * t42 + (-t11 - t24) * t39, t11 * t12 + t19 * t82 + t4 * t24 - t36 * t57; 0, 0, 0, 0, 0, 0, 0, t55, -0.2e1 * t61, t34, t25, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t69, -0.2e1 * t12 * t45 + 0.2e1 * t17, 0.2e1 * t82, -0.2e1 * t12 * t42 - 0.2e1 * t24 * t39, 0.2e1 * t24 * t12 + 0.2e1 * t36 * t82; 0, 0, 0, 0, 0, 0, 0, -t9, t8, t34, t25, 0, 0, 0, t56 - t65, -t64 + t73, t74 * t45 + t10 + t23, -t57, t74 * t42 + (-t11 - t29) * t39, -pkin(8) * t57 + t11 * t22 + t4 * t29; 0, 0, 0, 0, 0, 0, 0, -t62, -t61, t34, t25, 0, 0, 0, t49 - t65, -t64 + t69, t70 * t45 + t17 + t23, t82, t70 * t42 + (-t24 - t29) * t39, pkin(8) * t82 + t12 * t29 + t24 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t25, 0, 0, 0, -0.2e1 * t65, -0.2e1 * t64, -0.2e1 * t22 * t45 + 0.2e1 * t23, 0, -0.2e1 * t22 * t42 - 0.2e1 * t29 * t39, 0.2e1 * t29 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t66, 0, -t3, t2, -t3, t21, -t2, t21 * t19 + t51 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t66, 0, -t16, t15, -t16, t21, -t15, t21 * t36 - t51 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t66, 0, -t59, t60, -t59, t21, -t60, t21 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, qJ(5) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
