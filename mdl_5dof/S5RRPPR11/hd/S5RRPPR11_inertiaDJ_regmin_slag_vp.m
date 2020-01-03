% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:01
% DurationCPUTime: 0.61s
% Computational Cost: add. (517->110), mult. (1275->222), div. (0->0), fcn. (1024->6), ass. (0->73)
t89 = pkin(3) + pkin(6);
t56 = -pkin(2) - qJ(4);
t60 = cos(qJ(2));
t58 = sin(qJ(2));
t80 = t58 * qJ(3);
t64 = -t56 * t60 + t80;
t54 = sin(pkin(8));
t55 = cos(pkin(8));
t57 = sin(qJ(5));
t59 = cos(qJ(5));
t92 = -t57 * t54 + t59 * t55;
t38 = (t54 ^ 2 + t55 ^ 2) * qJD(4);
t77 = t60 * qJD(3);
t91 = t64 * qJD(2) + qJD(4) * t58 - t77;
t90 = 0.2e1 * qJD(3);
t88 = -pkin(7) + t56;
t78 = t58 * qJD(2);
t34 = t89 * t78;
t87 = t34 * t54;
t86 = t54 * t58;
t85 = t55 * t60;
t71 = pkin(2) * t78 - t58 * qJD(3);
t18 = -t60 * qJD(4) + (-qJ(3) * t60 + qJ(4) * t58) * qJD(2) + t71;
t49 = t60 * qJD(2);
t47 = pkin(6) * t49;
t35 = pkin(3) * t49 + t47;
t9 = t55 * t18 + t54 * t35;
t29 = -pkin(1) - t64;
t41 = t89 * t58;
t16 = t55 * t29 + t54 * t41;
t42 = t89 * t60;
t79 = qJD(5) * t60;
t76 = qJ(3) * qJD(3);
t75 = -0.2e1 * pkin(1) * qJD(2);
t74 = pkin(6) * t78;
t73 = t54 * t78;
t72 = t55 * t78;
t8 = -t54 * t18 + t55 * t35;
t3 = t9 * t54 + t8 * t55;
t70 = -t60 * pkin(2) - t80;
t33 = t55 * t41;
t10 = t58 * pkin(4) + t33 + (pkin(7) * t60 - t29) * t54;
t13 = -pkin(7) * t85 + t16;
t69 = t59 * t10 - t57 * t13;
t68 = t57 * t10 + t59 * t13;
t36 = t88 * t54;
t37 = t88 * t55;
t67 = t59 * t36 + t57 * t37;
t66 = t57 * t36 - t59 * t37;
t30 = t59 * t54 + t57 * t55;
t26 = t92 * qJD(5);
t62 = -t26 * t58 - t30 * t49;
t61 = t70 * qJD(2) + t77;
t46 = t54 * pkin(4) + qJ(3);
t44 = 0.2e1 * t58 * t49;
t39 = -pkin(1) + t70;
t25 = t30 * qJD(5);
t24 = pkin(4) * t85 + t42;
t23 = -qJ(3) * t49 + t71;
t21 = t30 * t60;
t20 = t92 * t60;
t19 = (-pkin(4) * t55 - t89) * t78;
t15 = -t54 * t29 + t33;
t14 = -t25 * t58 + t49 * t92;
t12 = t30 * t79 - t57 * t73 + t59 * t72;
t11 = t30 * t78 - t79 * t92;
t7 = -qJD(4) * t92 - t67 * qJD(5);
t6 = t30 * qJD(4) + t66 * qJD(5);
t5 = pkin(7) * t72 + t9;
t4 = (pkin(4) * t60 - pkin(7) * t86) * qJD(2) + t8;
t2 = -t68 * qJD(5) + t59 * t4 - t57 * t5;
t1 = -t69 * qJD(5) - t57 * t4 - t59 * t5;
t17 = [0, 0, 0, t44, 0.2e1 * (-t58 ^ 2 + t60 ^ 2) * qJD(2), 0, 0, 0, t58 * t75, t60 * t75, 0, 0.2e1 * t23 * t60 - 0.2e1 * t39 * t78, -0.2e1 * t23 * t58 - 0.2e1 * t39 * t49, 0.2e1 * t39 * t23, -0.2e1 * t34 * t85 + 0.2e1 * t8 * t58 + 0.2e1 * (-t42 * t55 * t58 + t15 * t60) * qJD(2), 0.2e1 * t60 * t87 - 0.2e1 * t9 * t58 + 0.2e1 * (-t16 * t60 + t42 * t86) * qJD(2), 0.2e1 * (t54 * t8 - t55 * t9) * t60 + 0.2e1 * (-t15 * t54 + t16 * t55) * t78, 0.2e1 * t15 * t8 + 0.2e1 * t16 * t9 - 0.2e1 * t42 * t34, -0.2e1 * t21 * t11, -0.2e1 * t11 * t20 - 0.2e1 * t21 * t12, 0.2e1 * t11 * t58 - 0.2e1 * t21 * t49, 0.2e1 * t12 * t58 - 0.2e1 * t20 * t49, t44, -0.2e1 * t24 * t12 + 0.2e1 * t19 * t20 + 0.2e1 * t2 * t58 + 0.2e1 * t69 * t49, 0.2e1 * t1 * t58 + 0.2e1 * t24 * t11 - 0.2e1 * t19 * t21 - 0.2e1 * t68 * t49; 0, 0, 0, 0, 0, t49, -t78, 0, -t47, t74, t61, t47, -t74, t61 * pkin(6), -t91 * t55 - t87, -t34 * t55 + t91 * t54, -t3, -t34 * qJ(3) + t42 * qJD(3) + t3 * t56 + (-t15 * t55 - t16 * t54) * qJD(4), t11 * t92 + t21 * t25, -t11 * t30 + t12 * t92 + t25 * t20 + t21 * t26, t14, t62, 0, qJD(3) * t20 - t46 * t12 + t19 * t30 + t24 * t26 - t66 * t49 + t7 * t58, -qJD(3) * t21 + t46 * t11 + t19 * t92 - t24 * t25 - t67 * t49 + t6 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0.2e1 * t76, t54 * t90, t55 * t90, 0.2e1 * t38, -0.2e1 * t56 * t38 + 0.2e1 * t76, -0.2e1 * t92 * t25, 0.2e1 * t25 * t30 - 0.2e1 * t26 * t92, 0, 0, 0, 0.2e1 * qJD(3) * t30 + 0.2e1 * t46 * t26, 0.2e1 * qJD(3) * t92 - 0.2e1 * t46 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, t47, t55 * t49, -t54 * t49, 0, t3, 0, 0, 0, 0, 0, t14, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t73, 0, -t34, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, t26, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, t49, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
