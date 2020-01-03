% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:05
% EndTime: 2019-12-31 20:11:07
% DurationCPUTime: 0.58s
% Computational Cost: add. (485->108), mult. (1073->203), div. (0->0), fcn. (741->4), ass. (0->78)
t87 = pkin(3) + pkin(6);
t46 = cos(qJ(2));
t47 = -pkin(2) - pkin(7);
t44 = sin(qJ(2));
t80 = t44 * qJ(3);
t53 = -t46 * t47 + t80;
t20 = -pkin(1) - t53;
t29 = t87 * t44;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t85 = t45 * t20 + t43 * t29;
t42 = t46 ^ 2;
t59 = qJD(2) * (t44 ^ 2 - t42);
t39 = t43 ^ 2;
t84 = -t45 ^ 2 + t39;
t58 = t84 * qJD(4);
t70 = t46 * qJD(3);
t30 = t87 * t46;
t73 = t30 * qJD(4);
t89 = t53 * qJD(2) - t70 - t73;
t88 = 0.2e1 * qJD(3);
t86 = pkin(4) * t45;
t82 = qJ(3) * t46;
t81 = qJ(5) * t46;
t79 = qJ(5) - t47;
t78 = qJD(2) * t30;
t77 = qJD(4) * t43;
t76 = qJD(4) * t45;
t75 = qJD(4) * t46;
t74 = qJD(4) * t47;
t72 = t44 * qJD(2);
t71 = t45 * qJD(5);
t36 = t46 * qJD(2);
t69 = qJ(3) * qJD(4);
t68 = -0.2e1 * pkin(1) * qJD(2);
t67 = pkin(4) * t77;
t66 = pkin(6) * t72;
t65 = t43 * t75;
t64 = t45 * t75;
t63 = t44 * t36;
t62 = t45 * t72;
t61 = t43 * t76;
t27 = t79 * t45;
t60 = -t20 + t81;
t57 = pkin(2) * t72 - t44 * qJD(3);
t56 = t43 * t62;
t23 = t45 * t29;
t6 = t44 * pkin(4) + t60 * t43 + t23;
t7 = -t45 * t81 + t85;
t55 = -t43 * t6 + t45 * t7;
t54 = -t46 * pkin(2) - t80;
t10 = (pkin(7) * t44 - t82) * qJD(2) + t57;
t34 = pkin(6) * t36;
t25 = pkin(3) * t36 + t34;
t4 = -t45 * t10 + t20 * t77 - t43 * t25 - t29 * t76;
t51 = t62 + t65;
t50 = t43 * t72 - t64;
t24 = t87 * t72;
t49 = -t24 + (-t44 * t47 - t82) * qJD(4);
t48 = t54 * qJD(2) + t70;
t33 = t43 * pkin(4) + qJ(3);
t32 = pkin(4) * t76 + qJD(3);
t31 = 0.2e1 * t63;
t28 = -pkin(1) + t54;
t26 = t79 * t43;
t19 = t45 * t25;
t17 = -t43 * t36 - t44 * t76;
t16 = t45 * t36 - t44 * t77;
t15 = t46 * t86 + t30;
t14 = -qJ(3) * t36 + t57;
t12 = -qJD(4) * t27 - t43 * qJD(5);
t11 = t79 * t77 - t71;
t8 = -pkin(4) * t65 + (-t86 - t87) * t72;
t5 = -t85 * qJD(4) - t43 * t10 + t19;
t3 = t11 * t45 + t12 * t43 + (-t26 * t45 + t27 * t43) * qJD(4);
t2 = t51 * qJ(5) - t46 * t71 - t4;
t1 = pkin(4) * t36 + t19 + t60 * t76 + (-qJ(5) * t72 - qJD(4) * t29 + qJD(5) * t46 - t10) * t43;
t9 = [0, 0, 0, t31, -0.2e1 * t59, 0, 0, 0, t44 * t68, t46 * t68, 0, 0.2e1 * t14 * t46 - 0.2e1 * t28 * t72, -0.2e1 * t14 * t44 - 0.2e1 * t28 * t36, 0.2e1 * t28 * t14, -0.2e1 * t39 * t63 + 0.2e1 * t42 * t61, -0.2e1 * t42 * t58 - 0.4e1 * t46 * t56, 0.2e1 * t43 * t59 - 0.2e1 * t44 * t64, 0.2e1 * t44 * t65 + 0.2e1 * t45 * t59, t31, 0.2e1 * (-t45 * t78 + t5) * t44 + 0.2e1 * ((-t43 * t20 + t23) * qJD(2) - t24 * t45 - t43 * t73) * t46, 0.2e1 * (t43 * t78 + t4) * t44 + 0.2e1 * (-t85 * qJD(2) + t24 * t43 - t45 * t73) * t46, 0.2e1 * t55 * t72 + 0.2e1 * (t1 * t43 - t2 * t45 + (t43 * t7 + t45 * t6) * qJD(4)) * t46, 0.2e1 * t6 * t1 + 0.2e1 * t15 * t8 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, t36, -t72, 0, -t34, t66, t48, t34, -t66, t48 * pkin(6), t46 * t58 + t56, 0.4e1 * t46 * t61 - t84 * t72, t16, t17, 0, t49 * t43 - t89 * t45, t89 * t43 + t49 * t45, (-t26 * t72 - t12 * t46 - t1 + (-t27 * t46 - t7) * qJD(4)) * t45 + (t27 * t72 + t11 * t46 - t2 + (-t26 * t46 + t6) * qJD(4)) * t43, -t1 * t27 + t6 * t11 + t7 * t12 + t15 * t32 - t2 * t26 + t8 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, qJ(3) * t88, -0.2e1 * t61, 0.2e1 * t58, 0, 0, 0, 0.2e1 * qJD(3) * t43 + 0.2e1 * t45 * t69, 0.2e1 * qJD(3) * t45 - 0.2e1 * t43 * t69, -0.2e1 * t3, -0.2e1 * t27 * t11 - 0.2e1 * t26 * t12 + 0.2e1 * t33 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, t34, 0, 0, 0, 0, 0, t16, t17, 0, t55 * qJD(4) + t1 * t45 + t2 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, t36, t5, t4, -t50 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, -t43 * t74, -t45 * t74, t67, t11 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
