% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:57
% EndTime: 2019-03-09 08:31:59
% DurationCPUTime: 0.74s
% Computational Cost: add. (1435->146), mult. (3045->256), div. (0->0), fcn. (2781->6), ass. (0->89)
t57 = sin(pkin(9));
t58 = cos(pkin(9));
t60 = sin(qJ(2));
t62 = cos(qJ(2));
t42 = t57 * t60 - t58 * t62;
t43 = t57 * t62 + t58 * t60;
t79 = -t62 * pkin(2) - pkin(1);
t72 = -t43 * qJ(4) + t79;
t99 = pkin(3) + pkin(8);
t18 = t99 * t42 + t72;
t94 = -qJ(3) - pkin(7);
t46 = t94 * t60;
t47 = t94 * t62;
t29 = -t58 * t46 - t57 * t47;
t24 = t43 * pkin(4) + t29;
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t93 = t61 * t18 + t59 * t24;
t55 = t59 ^ 2;
t56 = t61 ^ 2;
t76 = (t55 - t56) * qJD(5);
t100 = 2 * qJD(4);
t77 = qJD(2) * t94;
t35 = t62 * qJD(3) + t60 * t77;
t36 = -t60 * qJD(3) + t62 * t77;
t23 = t58 * t35 + t57 * t36;
t37 = t43 * qJD(2);
t14 = -t37 * pkin(4) + t23;
t98 = t14 * t42;
t97 = t55 * t37;
t96 = t59 * t37;
t95 = t61 * t37;
t91 = qJ(6) * t42;
t53 = -t58 * pkin(2) - pkin(3);
t50 = -pkin(8) + t53;
t90 = qJ(6) - t50;
t89 = qJD(5) * t59;
t88 = qJD(5) * t61;
t87 = t60 * qJD(2);
t86 = t61 * qJD(6);
t85 = t62 * qJD(2);
t84 = -0.2e1 * pkin(1) * qJD(2);
t83 = t59 * t95;
t54 = pkin(2) * t87;
t82 = pkin(5) * t89;
t81 = pkin(5) * t88;
t80 = t59 * t88;
t51 = t57 * pkin(2) + qJ(4);
t22 = t57 * t35 - t58 * t36;
t40 = t90 * t61;
t78 = -t18 - t91;
t21 = t61 * t24;
t5 = t43 * pkin(5) + t78 * t59 + t21;
t6 = t61 * t91 + t93;
t75 = -t5 * t59 + t6 * t61;
t30 = t57 * t46 - t58 * t47;
t74 = t29 * t22 + t30 * t23;
t73 = -qJD(4) * t42 - t51 * t37;
t71 = t42 * t88 + t96;
t70 = t42 * t89 - t95;
t38 = -t57 * t87 + t58 * t85;
t69 = t59 * t38 + t43 * t88;
t68 = -t61 * t38 + t43 * t89;
t67 = -t38 * qJ(4) - t43 * qJD(4) + t54;
t10 = t99 * t37 + t67;
t13 = t38 * pkin(4) + t22;
t3 = -t61 * t10 - t59 * t13 + t18 * t89 - t24 * t88;
t66 = t14 + (t42 * t51 - t43 * t50) * qJD(5);
t25 = -t42 * pkin(4) + t30;
t65 = -qJD(5) * t25 - t38 * t50 - t73;
t12 = t61 * t13;
t1 = t38 * pkin(5) + t12 + t78 * t88 + (-qJ(6) * t37 - qJD(5) * t24 - qJD(6) * t42 - t10) * t59;
t2 = -t70 * qJ(6) + t42 * t86 - t3;
t64 = -t1 * t59 + t2 * t61 + (-t5 * t61 - t59 * t6) * qJD(5);
t63 = 0.2e1 * t22 * t43 - 0.2e1 * t23 * t42 + 0.2e1 * t29 * t38 - 0.2e1 * t30 * t37;
t48 = qJD(4) + t81;
t45 = t59 * pkin(5) + t51;
t41 = t42 ^ 2;
t39 = t90 * t59;
t34 = t56 * t37;
t32 = -qJD(5) * t40 - t59 * qJD(6);
t31 = t90 * t89 - t86;
t28 = t42 * pkin(3) + t72;
t16 = t37 * pkin(3) + t67;
t15 = (-pkin(5) * t61 - pkin(4)) * t42 + t30;
t8 = t31 * t61 + t32 * t59 + (-t39 * t61 + t40 * t59) * qJD(5);
t7 = t70 * pkin(5) + t14;
t4 = -t93 * qJD(5) - t59 * t10 + t12;
t9 = [0, 0, 0, 0.2e1 * t60 * t85, 0.2e1 * (-t60 ^ 2 + t62 ^ 2) * qJD(2), 0, 0, 0, t60 * t84, t62 * t84, t63, 0.2e1 * t79 * t54 + 0.2e1 * t74, t63, -0.2e1 * t16 * t42 - 0.2e1 * t28 * t37, -0.2e1 * t16 * t43 - 0.2e1 * t28 * t38, 0.2e1 * t28 * t16 + 0.2e1 * t74, 0.2e1 * t41 * t80 + 0.2e1 * t42 * t97, -0.2e1 * t41 * t76 + 0.4e1 * t42 * t83, 0.2e1 * t69 * t42 + 0.2e1 * t43 * t96, -0.2e1 * t68 * t42 + 0.2e1 * t43 * t95, 0.2e1 * t43 * t38, 0.2e1 * t4 * t43 + 0.2e1 * (-t59 * t18 + t21) * t38 - 0.2e1 * t61 * t98 + 0.2e1 * t70 * t25, 0.2e1 * t71 * t25 + 0.2e1 * t3 * t43 - 0.2e1 * t93 * t38 + 0.2e1 * t59 * t98, 0.2e1 * t75 * t37 + 0.2e1 * t64 * t42, 0.2e1 * t5 * t1 + 0.2e1 * t15 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t85, -t87, 0, -pkin(7) * t85, pkin(7) * t87 (-t37 * t57 - t38 * t58) * pkin(2) (-t22 * t58 + t23 * t57) * pkin(2), t53 * t38 + t73, t22, t23, t30 * qJD(4) + t22 * t53 + t23 * t51, -t42 * t76 + t83, -0.4e1 * t42 * t80 + t34 - t97, -t68, -t69, 0, t66 * t59 - t65 * t61, t65 * t59 + t66 * t61 (t32 * t42 - t37 * t39 - t1 + (t40 * t42 - t6) * qJD(5)) * t61 + (-t31 * t42 + t37 * t40 - t2 + (t39 * t42 + t5) * qJD(5)) * t59, -t1 * t40 + t15 * t48 - t2 * t39 + t5 * t31 + t6 * t32 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t51 * t100, -0.2e1 * t80, 0.2e1 * t76, 0, 0, 0, 0.2e1 * qJD(4) * t59 + 0.2e1 * t51 * t88, 0.2e1 * qJD(4) * t61 - 0.2e1 * t51 * t89, -0.2e1 * t8, -0.2e1 * t40 * t31 - 0.2e1 * t39 * t32 + 0.2e1 * t45 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t37, -t38, t16, 0, 0, 0, 0, 0, -t69, t68, t34 + t97, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t59 + t32 * t61 + (t39 * t59 + t40 * t61) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, t22, 0, 0, 0, 0, 0, -t68, -t69, 0, t75 * qJD(5) + t1 * t61 + t2 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t70, t38, t4, t3, -t71 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0, -t50 * t89, -t50 * t88, t82, t31 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t89, 0, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
