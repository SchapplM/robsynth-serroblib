% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:54
% EndTime: 2019-03-09 02:13:56
% DurationCPUTime: 0.74s
% Computational Cost: add. (1192->132), mult. (2433->241), div. (0->0), fcn. (2285->6), ass. (0->83)
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t98 = sin(qJ(4));
t99 = cos(qJ(4));
t33 = t99 * t51 + t98 * t52;
t43 = t51 * pkin(3) + qJ(2);
t58 = -t98 * t51 + t99 * t52;
t20 = t33 * pkin(4) - pkin(8) * t58 + t43;
t53 = -pkin(1) - qJ(3);
t100 = -pkin(7) + t53;
t35 = t100 * t51;
t36 = t100 * t52;
t23 = t99 * t35 + t98 * t36;
t55 = cos(qJ(5));
t21 = t55 * t23;
t54 = sin(qJ(5));
t90 = t54 * t20 + t21;
t49 = t54 ^ 2;
t50 = t55 ^ 2;
t88 = t49 - t50;
t68 = qJD(5) * t88;
t37 = (t51 ^ 2 + t52 ^ 2) * qJD(3);
t31 = t58 ^ 2;
t101 = 2 * qJD(2);
t70 = qJD(4) * t98;
t71 = qJD(4) * t99;
t28 = -t51 * t71 - t52 * t70;
t97 = t58 * t28;
t96 = t58 * t54;
t95 = t58 * t55;
t29 = -t51 * t70 + t52 * t71;
t94 = t33 * t29;
t93 = t54 * t28;
t92 = t55 * t28;
t91 = -qJ(6) - pkin(8);
t87 = t49 + t50;
t86 = qJ(6) * t58;
t85 = qJD(5) * t54;
t84 = qJD(5) * t55;
t83 = qJ(2) * qJD(2);
t82 = -0.2e1 * pkin(4) * qJD(5);
t19 = t29 * pkin(4) - t28 * pkin(8) + qJD(2);
t9 = t33 * qJD(3) + t35 * t70 - t36 * t71;
t81 = t54 * t19 + t20 * t84 - t55 * t9;
t80 = pkin(5) * t85;
t79 = t58 * t85;
t78 = t58 * t84;
t77 = t54 * t84;
t76 = t33 ^ 2 + t31;
t75 = t55 * t19 + t54 * t9;
t74 = t87 * t29;
t73 = -0.4e1 * t54 * t95;
t72 = t55 * t20 - t54 * t23;
t69 = qJD(5) * t91;
t67 = -pkin(4) * t28 - pkin(8) * t29;
t66 = -pkin(4) * t58 - pkin(8) * t33;
t5 = t33 * pkin(5) - t55 * t86 + t72;
t6 = -t54 * t86 + t90;
t65 = -t5 * t55 - t54 * t6;
t64 = -t5 * t54 + t55 * t6;
t22 = t98 * t35 - t99 * t36;
t38 = t91 * t54;
t39 = t91 * t55;
t63 = -t38 * t54 - t39 * t55;
t62 = -qJ(6) * t28 - qJD(6) * t58;
t61 = -t78 - t93;
t60 = -t79 + t92;
t18 = t54 * t29 + t33 * t84;
t17 = -t55 * t29 + t33 * t85;
t59 = -0.2e1 * t94 - 0.2e1 * t97;
t1 = t29 * pkin(5) + t62 * t55 + (-t21 + (-t20 + t86) * t54) * qJD(5) + t75;
t2 = -qJ(6) * t78 + (-qJD(5) * t23 + t62) * t54 + t81;
t57 = t64 * qJD(5) + t1 * t55 + t2 * t54;
t25 = t55 * qJD(6) + t54 * t69;
t26 = -t54 * qJD(6) + t55 * t69;
t56 = t25 * t55 - t26 * t54 + (-t38 * t55 + t39 * t54) * qJD(5);
t10 = qJD(3) * t58 + qJD(4) * t23;
t45 = -t55 * pkin(5) - pkin(4);
t11 = pkin(5) * t96 + t22;
t7 = -t61 * pkin(5) + t10;
t4 = -qJD(5) * t90 + t75;
t3 = t23 * t85 - t81;
t8 = [0, 0, 0, 0, t101, 0.2e1 * t83, t51 * t101, t52 * t101, 0.2e1 * t37, -0.2e1 * t37 * t53 + 0.2e1 * t83, 0.2e1 * t97, -0.2e1 * t28 * t33 - 0.2e1 * t29 * t58, 0, 0, 0, 0.2e1 * qJD(2) * t33 + 0.2e1 * t43 * t29, 0.2e1 * qJD(2) * t58 + 0.2e1 * t43 * t28, -0.2e1 * t31 * t77 + 0.2e1 * t50 * t97, t28 * t73 + 0.2e1 * t31 * t68, -0.2e1 * t17 * t58 + 0.2e1 * t33 * t92, -0.2e1 * t18 * t58 - 0.2e1 * t33 * t93, 0.2e1 * t94, 0.2e1 * t10 * t96 - 0.2e1 * t61 * t22 + 0.2e1 * t72 * t29 + 0.2e1 * t4 * t33, 0.2e1 * t10 * t95 + 0.2e1 * t60 * t22 - 0.2e1 * t90 * t29 + 0.2e1 * t3 * t33, 0.2e1 * t65 * t28 - 0.2e1 * t57 * t58, 0.2e1 * t5 * t1 + 0.2e1 * t11 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t59 - t76 * t84, t55 * t59 + t76 * t85, 0, -t11 * t28 - t7 * t58 + t64 * t29 + (t65 * qJD(5) - t1 * t54 + t2 * t55) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t33 * t74 + 0.2e1 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t29, t28, 0, 0, 0, 0, 0, -t17, -t18, -t87 * t28, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, -t10, t9, t54 * t92 - t58 * t68, qJD(5) * t73 - t88 * t28, t18, -t17, 0, -t10 * t55 + t67 * t54 + (t22 * t54 + t66 * t55) * qJD(5), t10 * t54 + t67 * t55 + (t22 * t55 - t66 * t54) * qJD(5) (-t26 * t58 - t28 * t38 + t2 + (t39 * t58 - t5) * qJD(5)) * t55 + (-t25 * t58 + t28 * t39 - t1 + (t38 * t58 - t6) * qJD(5)) * t54, t1 * t38 + t11 * t80 - t2 * t39 + t6 * t25 + t5 * t26 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, 0, 0, 0, 0, t60, t61, t74, -pkin(5) * t79 - t28 * t45 + t63 * t29 + t56 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * qJD(5) + t54 * t25 + t55 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t77, -0.2e1 * t68, 0, 0, 0, t54 * t82, t55 * t82, 0.2e1 * t56, -0.2e1 * t39 * t25 + 0.2e1 * t38 * t26 + 0.2e1 * t45 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t61, t29, t4, t3, -t60 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -t18 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t84, 0, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t85, 0, -pkin(8) * t84, pkin(8) * t85, -pkin(5) * t84, t26 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
