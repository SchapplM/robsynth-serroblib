% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP2
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
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:58
% EndTime: 2019-12-31 20:52:00
% DurationCPUTime: 0.37s
% Computational Cost: add. (366->96), mult. (831->143), div. (0->0), fcn. (503->4), ass. (0->64)
t53 = sin(qJ(3));
t47 = t53 * qJ(4);
t55 = cos(qJ(3));
t79 = t55 * pkin(3) + t47;
t58 = 2 * qJD(4);
t56 = cos(qJ(2));
t42 = -pkin(1) * t56 - pkin(2);
t24 = t42 - t79;
t48 = t55 * pkin(4);
t16 = -t24 + t48;
t46 = t55 * qJD(3);
t54 = sin(qJ(2));
t72 = pkin(1) * qJD(2);
t64 = t54 * t72;
t70 = t53 * qJD(3);
t22 = pkin(3) * t70 - qJ(4) * t46 - t53 * qJD(4);
t9 = -pkin(4) * t70 - t22;
t5 = t9 - t64;
t78 = t16 * t46 + t5 * t53;
t67 = pkin(2) + t79;
t25 = t48 + t67;
t77 = t25 * t46 + t9 * t53;
t76 = pkin(7) - qJ(5);
t10 = t64 + t22;
t75 = -t10 - t22;
t74 = t42 * t46 + t53 * t64;
t68 = qJ(5) * qJD(3);
t73 = -t53 * qJD(5) - t55 * t68;
t41 = pkin(1) * t54 + pkin(7);
t71 = -qJ(5) + t41;
t69 = t55 * qJD(4);
t66 = pkin(2) * t70;
t65 = pkin(2) * t46;
t63 = t56 * t72;
t62 = pkin(7) * t70;
t61 = pkin(7) * t46;
t60 = -qJD(5) * t55 + t53 * t68;
t11 = t41 * t70 - t55 * t63;
t51 = t53 ^ 2;
t52 = t55 ^ 2;
t18 = (t51 + t52) * t63;
t59 = t42 * t70 - t55 * t64;
t12 = t41 * t46 + t53 * t63;
t21 = -qJD(3) * t79 + t69;
t57 = -pkin(3) - pkin(4);
t50 = qJ(4) * t58;
t37 = 0.2e1 * t53 * t46;
t34 = t76 * t55;
t33 = t76 * t53;
t29 = 0.2e1 * (-t51 + t52) * qJD(3);
t28 = t34 * t70;
t27 = t71 * t55;
t26 = t71 * t53;
t23 = t67 * t70;
t20 = t61 + t73;
t19 = t60 - t62;
t17 = -t69 + (-t55 * t57 + t47) * qJD(3);
t15 = t27 * t70;
t13 = t24 * t70;
t7 = t9 * t55;
t4 = t12 + t73;
t3 = -t11 + t60;
t2 = t5 * t55;
t1 = [0, 0, 0, 0, -0.2e1 * t64, -0.2e1 * t63, t37, t29, 0, 0, 0, 0.2e1 * t59, 0.2e1 * t74, -0.2e1 * t10 * t55 + 0.2e1 * t13, 0.2e1 * t18, -0.2e1 * t10 * t53 - 0.2e1 * t24 * t46, 0.2e1 * t24 * t10 + 0.2e1 * t18 * t41, -0.2e1 * t16 * t70 + 0.2e1 * t2, 0.2e1 * t78, -0.2e1 * t4 * t53 + 0.2e1 * t15 + 0.2e1 * (-qJD(3) * t26 - t3) * t55, 0.2e1 * t16 * t5 + 0.2e1 * t26 * t4 + 0.2e1 * t27 * t3; 0, 0, 0, 0, -t64, -t63, t37, t29, 0, 0, 0, t59 - t66, -t65 + t74, t55 * t75 + t13 - t23, t18, t75 * t53 + (-t24 + t67) * t46, pkin(7) * t18 - t10 * t67 + t24 * t22, t2 + t7 + (-t16 - t25) * t70, t77 + t78, t15 + t28 + (-t20 - t4) * t53 + (-t19 - t3 + (-t26 - t33) * qJD(3)) * t55, t16 * t9 + t19 * t27 + t20 * t26 + t25 * t5 + t3 * t34 + t33 * t4; 0, 0, 0, 0, 0, 0, t37, t29, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t65, -0.2e1 * t22 * t55 - 0.2e1 * t23, 0, -0.2e1 * t22 * t53 + 0.2e1 * t46 * t67, -0.2e1 * t67 * t22, -0.2e1 * t25 * t70 + 0.2e1 * t7, 0.2e1 * t77, -0.2e1 * t20 * t53 + 0.2e1 * t28 + 0.2e1 * (-qJD(3) * t33 - t19) * t55, 0.2e1 * t19 * t34 + 0.2e1 * t20 * t33 + 0.2e1 * t25 * t9; 0, 0, 0, 0, 0, 0, 0, 0, t46, -t70, 0, -t12, t11, -t12, t21, -t11, (-pkin(3) * t53 + qJ(4) * t55) * t63 + t21 * t41, -t4, t3, t17, qJ(4) * t3 + qJD(4) * t27 + t4 * t57; 0, 0, 0, 0, 0, 0, 0, 0, t46, -t70, 0, -t61, t62, -t61, t21, -t62, t21 * pkin(7), -t20, t19, t17, qJ(4) * t19 + qJD(4) * t34 + t20 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t50, 0, t58, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t12, 0, 0, -t46, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t61, 0, 0, -t46, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t46, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t46, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
