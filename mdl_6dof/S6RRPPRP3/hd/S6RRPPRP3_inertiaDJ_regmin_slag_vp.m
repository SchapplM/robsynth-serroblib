% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:33
% EndTime: 2019-03-09 08:35:35
% DurationCPUTime: 0.72s
% Computational Cost: add. (770->138), mult. (1488->248), div. (0->0), fcn. (1036->4), ass. (0->86)
t54 = sin(qJ(2));
t42 = t54 * qJ(3);
t56 = cos(qJ(2));
t98 = -pkin(2) * t56 - t42;
t30 = -pkin(1) + t98;
t26 = pkin(3) * t56 - t30;
t15 = pkin(4) * t54 + pkin(8) * t56 + t26;
t31 = (pkin(7) - qJ(4)) * t54;
t55 = cos(qJ(5));
t27 = t55 * t31;
t53 = sin(qJ(5));
t95 = t15 * t53 + t27;
t48 = t53 ^ 2;
t50 = t55 ^ 2;
t93 = t48 - t50;
t69 = t93 * qJD(5);
t85 = t54 * qJD(2);
t77 = pkin(7) * t85;
t82 = qJ(4) * qJD(2);
t20 = t56 * qJD(4) - t54 * t82 + t77;
t52 = qJ(3) + pkin(4);
t57 = -pkin(2) - pkin(3);
t47 = -pkin(8) + t57;
t96 = t47 * t54;
t97 = (t52 * t56 + t96) * qJD(5) + t20;
t58 = 2 * qJD(3);
t40 = t56 * qJD(2);
t94 = qJ(3) * t40 + qJD(3) * t54;
t51 = t56 ^ 2;
t92 = t54 ^ 2 - t51;
t91 = qJ(6) * t56;
t90 = t55 * qJ(6);
t89 = qJ(6) - t47;
t88 = qJD(2) * t55;
t41 = qJD(5) * t53;
t87 = qJD(5) * t55;
t86 = qJD(5) * t56;
t84 = t55 * qJD(6);
t83 = t56 * qJD(3);
t81 = -0.2e1 * pkin(1) * qJD(2);
t10 = (pkin(4) * t56 + t96) * qJD(2) + t94;
t38 = pkin(7) * t40;
t21 = -qJD(4) * t54 - t56 * t82 + t38;
t80 = t10 * t53 + t15 * t87 + t21 * t55;
t79 = pkin(5) * t41;
t78 = pkin(5) * t87;
t76 = t53 * t86;
t75 = t55 * t86;
t74 = t53 * t87;
t73 = t54 * t40;
t72 = t55 * t85;
t71 = t10 * t55 - t53 * t21;
t70 = t15 * t55 - t31 * t53;
t68 = t92 * qJD(2);
t67 = qJD(5) * t89;
t66 = t53 * t72;
t6 = pkin(5) * t54 + t56 * t90 + t70;
t7 = t53 * t91 + t95;
t65 = -t53 * t7 - t55 * t6;
t63 = t72 + t76;
t62 = t53 * t85 - t75;
t61 = qJD(2) * t98 + t83;
t1 = t56 * t84 + (pkin(5) * t56 - t54 * t90) * qJD(2) + (-t27 + (-t15 - t91) * t53) * qJD(5) + t71;
t2 = qJ(6) * t75 + (-qJ(6) * t85 - qJD(5) * t31 + qJD(6) * t56) * t53 + t80;
t60 = t1 * t55 + t2 * t53 + (-t53 * t6 + t55 * t7) * qJD(5);
t43 = t56 * qJ(4);
t32 = pkin(7) * t56 - t43;
t59 = -t83 - qJD(5) * t32 + (-t47 * t56 + t52 * t54) * qJD(2);
t46 = qJ(3) * t58;
t35 = pkin(5) * t55 + t52;
t34 = qJD(3) - t79;
t33 = 0.2e1 * t73;
t29 = t89 * t55;
t28 = t89 * t53;
t25 = -t40 * t53 - t54 * t87;
t24 = -t40 * t55 + t41 * t54;
t23 = -t43 + (-pkin(5) * t53 + pkin(7)) * t56;
t22 = pkin(2) * t85 - t94;
t19 = t53 * qJD(6) + t55 * t67;
t18 = t53 * t67 - t84;
t16 = t57 * t85 + t94;
t11 = pkin(5) * t62 - t20;
t5 = t18 * t55 - t19 * t53 + (-t28 * t55 + t29 * t53) * qJD(5);
t4 = -qJD(5) * t95 + t71;
t3 = t31 * t41 - t80;
t8 = [0, 0, 0, t33, -0.2e1 * t68, 0, 0, 0, t54 * t81, t56 * t81, -0.2e1 * t22 * t56 + 0.2e1 * t30 * t85, 0, -0.2e1 * t22 * t54 - 0.2e1 * t30 * t40, 0.2e1 * t30 * t22, 0.2e1 * t16 * t54 + 0.2e1 * t26 * t40, -0.2e1 * t16 * t56 + 0.2e1 * t26 * t85, 0.2e1 * t20 * t56 - 0.2e1 * t21 * t54 + 0.2e1 * (-t31 * t56 + t32 * t54) * qJD(2), 0.2e1 * t16 * t26 - 0.2e1 * t20 * t32 + 0.2e1 * t21 * t31, -0.2e1 * t50 * t73 - 0.2e1 * t51 * t74, 0.2e1 * t51 * t69 + 0.4e1 * t56 * t66, 0.2e1 * t54 * t76 + 0.2e1 * t88 * t92, -0.2e1 * t53 * t68 + 0.2e1 * t54 * t75, t33, 0.2e1 * (qJD(2) * t32 * t53 + t4) * t54 + 0.2e1 * (qJD(2) * t70 + t20 * t53 - t32 * t87) * t56, 0.2e1 * (t32 * t88 + t3) * t54 + 0.2e1 * (-qJD(2) * t95 + t20 * t55 + t32 * t41) * t56, 0.2e1 * t56 * t60 + 0.2e1 * t65 * t85, 0.2e1 * t1 * t6 + 0.2e1 * t11 * t23 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, t40, -t85, 0, -t38, t77, -t38, t61, -t77, t61 * pkin(7), -t20, t21, -t83 + (-t56 * t57 + t42) * qJD(2), -qJ(3) * t20 + qJD(3) * t32 + t21 * t57, -t56 * t69 - t66, -0.4e1 * t56 * t74 + t85 * t93, t25, t24, 0, t53 * t59 - t55 * t97, t53 * t97 + t55 * t59 (-t28 * t85 + t19 * t56 - t2 + (-t29 * t56 + t6) * qJD(5)) * t55 + (t29 * t85 + t18 * t56 + t1 + (-t28 * t56 + t7) * qJD(5)) * t53, t1 * t28 + t11 * t35 + t18 * t7 + t19 * t6 - t2 * t29 + t23 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t46, t58, 0, 0, t46, 0.2e1 * t74, -0.2e1 * t69, 0, 0, 0, 0.2e1 * qJD(3) * t55 - 0.2e1 * t41 * t52, -0.2e1 * qJD(3) * t53 - 0.2e1 * t52 * t87, -0.2e1 * t5, -0.2e1 * t18 * t29 + 0.2e1 * t19 * t28 + 0.2e1 * t34 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t38, 0, 0, -t40, t21, 0, 0, 0, 0, 0, t25, t24, 0, qJD(5) * t65 - t1 * t53 + t2 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t85, 0, t16, 0, 0, 0, 0, 0, -t24, t25 (-t48 - t50) * t85, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t53 + t19 * t55 + (-t28 * t53 - t29 * t55) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t62, t40, t4, t3, -t63 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t41, 0, -t47 * t87, t47 * t41, t78, t19 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t41, 0, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t87, 0, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
