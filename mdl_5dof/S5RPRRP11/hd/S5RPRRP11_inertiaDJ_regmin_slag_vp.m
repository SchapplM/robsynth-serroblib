% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:35
% DurationCPUTime: 0.68s
% Computational Cost: add. (1070->110), mult. (2456->216), div. (0->0), fcn. (2278->6), ass. (0->77)
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t47 = sin(qJ(3));
t89 = cos(qJ(3));
t31 = t89 * t44 + t47 * t45;
t37 = -t45 * pkin(2) - pkin(1);
t68 = t89 * t45;
t53 = -t47 * t44 + t68;
t20 = -pkin(3) * t53 - t31 * pkin(7) + t37;
t81 = pkin(6) + qJ(2);
t33 = t81 * t44;
t34 = t81 * t45;
t23 = -t47 * t33 + t89 * t34;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t97 = t46 * t20 + t48 * t23;
t42 = t46 ^ 2;
t43 = t48 ^ 2;
t79 = t42 - t43;
t65 = t79 * qJD(4);
t69 = t89 * t33;
t10 = qJD(3) * t69 - qJD(2) * t68 + (qJD(2) * t44 + qJD(3) * t34) * t47;
t25 = t53 * qJD(3);
t26 = t31 * qJD(3);
t19 = t26 * pkin(3) - t25 * pkin(7);
t4 = -qJD(4) * t97 + t46 * t10 + t48 * t19;
t58 = t48 * pkin(4) + t46 * qJ(5);
t96 = qJD(4) * t58 - t48 * qJD(5);
t39 = qJD(4) * t48;
t77 = qJD(4) * t46;
t3 = t48 * t10 - t46 * t19 - t20 * t39 + t23 * t77;
t76 = t53 * qJD(5);
t78 = t26 * qJ(5);
t1 = -t3 - t76 + t78;
t90 = t26 * pkin(4);
t2 = -t4 - t90;
t6 = -qJ(5) * t53 + t97;
t56 = t48 * t20 - t46 * t23;
t7 = pkin(4) * t53 - t56;
t95 = t1 * t46 - t2 * t48 + (t46 * t7 + t48 * t6) * qJD(4);
t94 = 0.2e1 * t37;
t93 = 0.2e1 * qJD(5);
t92 = pkin(7) * t26;
t91 = pkin(7) * t53;
t11 = qJD(2) * t31 + qJD(3) * t23;
t88 = t11 * t46;
t87 = t31 * t25;
t86 = t31 * t48;
t85 = t46 * t25;
t84 = t46 * t26;
t83 = t48 * t25;
t82 = t48 * t26;
t75 = t46 * qJD(5);
t73 = -0.2e1 * pkin(3) * qJD(4);
t72 = pkin(7) * t77;
t71 = pkin(7) * t39;
t70 = t46 * t39;
t67 = -0.4e1 * t46 * t86;
t64 = 0.2e1 * (t44 ^ 2 + t45 ^ 2) * qJD(2);
t63 = -pkin(3) * t25 - t92;
t62 = pkin(3) * t31 - t91;
t59 = -t46 * t6 + t48 * t7;
t22 = t47 * t34 + t69;
t57 = pkin(4) * t46 - qJ(5) * t48;
t54 = t31 * t39 + t85;
t17 = -t31 * t77 + t83;
t18 = -t39 * t53 + t84;
t16 = t53 * t77 + t82;
t32 = -pkin(3) - t58;
t5 = t57 * t25 + t96 * t31 + t11;
t52 = -t5 + (t31 * t32 + t91) * qJD(4);
t24 = -pkin(4) * t77 + qJ(5) * t39 + t75;
t8 = t31 * t57 + t22;
t51 = -qJD(4) * t8 + t24 * t31 - t25 * t32 + t92;
t49 = qJD(4) * t59 + t1 * t48 + t2 * t46;
t28 = t31 ^ 2;
t9 = [0, 0, 0, 0, 0, t64, qJ(2) * t64, 0.2e1 * t87, 0.2e1 * t25 * t53 - 0.2e1 * t31 * t26, 0, 0, 0, t26 * t94, t25 * t94, -0.2e1 * t28 * t70 + 0.2e1 * t43 * t87, t25 * t67 + 0.2e1 * t28 * t65, -0.2e1 * t17 * t53 + 0.2e1 * t31 * t82, -0.2e1 * t31 * t84 + 0.2e1 * t53 * t54, -0.2e1 * t53 * t26, 0.2e1 * t22 * t54 + 0.2e1 * t26 * t56 + 0.2e1 * t31 * t88 - 0.2e1 * t4 * t53, 0.2e1 * t11 * t86 + 0.2e1 * t17 * t22 - 0.2e1 * t26 * t97 - 0.2e1 * t3 * t53, 0.2e1 * t8 * t85 + 0.2e1 * t2 * t53 - 0.2e1 * t7 * t26 + 0.2e1 * (t8 * t39 + t5 * t46) * t31, 0.2e1 * t59 * t25 - 0.2e1 * t95 * t31, -0.2e1 * t8 * t83 - 0.2e1 * t1 * t53 + 0.2e1 * t6 * t26 + 0.2e1 * (-t5 * t48 + t8 * t77) * t31, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t25, 0, 0, 0, 0, 0, t16, -t18, t16, (-t42 - t43) * t25, t18, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t26, 0, -t11, t10, -t31 * t65 + t46 * t83, qJD(4) * t67 - t79 * t25, t18, t16, 0, -t11 * t48 + t63 * t46 + (t22 * t46 - t48 * t62) * qJD(4), t88 + t63 * t48 + (t22 * t48 + t46 * t62) * qJD(4), -t46 * t51 + t48 * t52, t49, t46 * t52 + t48 * t51, pkin(7) * t49 - t8 * t24 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t70, -0.2e1 * t65, 0, 0, 0, t46 * t73, t48 * t73, 0.2e1 * t24 * t48 + 0.2e1 * t32 * t77, 0, 0.2e1 * t24 * t46 - 0.2e1 * t32 * t39, -0.2e1 * t32 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t54, t26, t4, t3, t4 + 0.2e1 * t90, -t58 * t25 + (qJD(4) * t57 - t75) * t31, -t3 - 0.2e1 * t76 + 0.2e1 * t78, -t2 * pkin(4) + t1 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t39, -t77, 0, t39, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t77, 0, -t71, t72, -t71, -t96, -t72, -t96 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, qJ(5) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
