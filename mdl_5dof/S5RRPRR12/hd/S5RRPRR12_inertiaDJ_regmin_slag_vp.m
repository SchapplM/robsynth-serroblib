% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:40
% DurationCPUTime: 0.68s
% Computational Cost: add. (639->115), mult. (1517->220), div. (0->0), fcn. (1290->6), ass. (0->85)
t51 = cos(qJ(5));
t47 = t51 ^ 2;
t48 = sin(qJ(5));
t83 = t48 ^ 2 - t47;
t99 = qJD(5) * t83;
t102 = 0.2e1 * t99;
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t101 = -t53 * pkin(2) - t50 * qJ(3);
t42 = t53 * qJD(2);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t77 = t50 * qJD(2);
t100 = -t49 * t42 + t52 * t77;
t88 = t50 * t49;
t27 = t53 * t52 + t88;
t28 = -t53 * t49 + t50 * t52;
t96 = -pkin(2) - pkin(3);
t71 = t52 * t96;
t29 = t49 * qJ(3) + pkin(4) - t71;
t57 = t52 * qJ(3) + t49 * t96;
t30 = -pkin(8) + t57;
t95 = pkin(6) - pkin(7);
t33 = t95 * t50;
t34 = t95 * t53;
t85 = t52 * t34;
t14 = t49 * t33 + t85;
t7 = t14 * qJD(4) + (-t95 * t88 - t85) * qJD(2);
t98 = -t7 + (t27 * t30 - t28 * t29) * qJD(5);
t97 = 2 * qJD(3);
t94 = t7 * t48;
t93 = t7 * t51;
t17 = t49 * qJD(3) + t57 * qJD(4);
t92 = t17 * t48;
t91 = t17 * t51;
t12 = (qJD(2) - qJD(4)) * t27;
t90 = t28 * t12;
t81 = qJD(4) * t52;
t82 = qJD(4) * t49;
t11 = t50 * t81 - t53 * t82 - t100;
t89 = t48 * t11;
t87 = t51 * t11;
t86 = t51 * t12;
t84 = qJ(3) * t42 + t50 * qJD(3);
t80 = qJD(5) * t48;
t79 = qJD(5) * t51;
t78 = qJD(5) * t52;
t76 = -0.2e1 * pkin(1) * qJD(2);
t75 = -0.2e1 * pkin(4) * qJD(5);
t31 = -pkin(1) + t101;
t74 = t48 * t86;
t73 = pkin(6) * t77;
t72 = pkin(6) * t42;
t68 = t48 * t79;
t67 = qJD(5) * (pkin(4) + t29);
t23 = t53 * pkin(3) - t31;
t65 = -pkin(4) * t12 - pkin(8) * t11;
t64 = pkin(4) * t28 + pkin(8) * t27;
t10 = t27 * pkin(4) - t28 * pkin(8) + t23;
t63 = t51 * t10 - t48 * t14;
t62 = t48 * t10 + t51 * t14;
t60 = t27 * t49 + t28 * t52;
t59 = t48 * t12 + t28 * t79;
t58 = -t28 * t80 + t86;
t15 = t96 * t77 + t84;
t56 = t101 * qJD(2) + t53 * qJD(3);
t55 = -t11 * t49 - t12 * t52 + (-t27 * t52 + t28 * t49) * qJD(4);
t13 = -t52 * t33 + t49 * t34;
t16 = qJ(3) * t82 - t52 * qJD(3) - qJD(4) * t71;
t54 = -qJD(5) * t13 - t11 * t30 + t12 * t29 + t16 * t27 + t17 * t28;
t35 = 0.2e1 * t68;
t26 = -0.2e1 * t99;
t25 = t28 ^ 2;
t20 = t48 * t78 + t51 * t82;
t19 = t48 * t82 - t51 * t78;
t18 = pkin(2) * t77 - t84;
t9 = t27 * t79 + t89;
t8 = t27 * t80 - t87;
t6 = t100 * t95 - t33 * t81 + t34 * t82;
t5 = t28 * t99 - t74;
t4 = t11 * pkin(4) - t12 * pkin(8) + t15;
t3 = t83 * t12 + 0.4e1 * t28 * t68;
t2 = -t62 * qJD(5) + t51 * t4 + t48 * t6;
t1 = -t63 * qJD(5) - t48 * t4 + t51 * t6;
t21 = [0, 0, 0, 0.2e1 * t50 * t42, 0.2e1 * (-t50 ^ 2 + t53 ^ 2) * qJD(2), 0, 0, 0, t50 * t76, t53 * t76, -0.2e1 * t18 * t53 + 0.2e1 * t31 * t77, 0, -0.2e1 * t18 * t50 - 0.2e1 * t31 * t42, 0.2e1 * t31 * t18, 0.2e1 * t90, -0.2e1 * t28 * t11 - 0.2e1 * t12 * t27, 0, 0, 0, 0.2e1 * t23 * t11 + 0.2e1 * t15 * t27, 0.2e1 * t23 * t12 + 0.2e1 * t15 * t28, -0.2e1 * t25 * t68 + 0.2e1 * t47 * t90, t25 * t102 - 0.4e1 * t28 * t74, 0.2e1 * t58 * t27 + 0.2e1 * t28 * t87, -0.2e1 * t59 * t27 - 0.2e1 * t28 * t89, 0.2e1 * t27 * t11, 0.2e1 * t63 * t11 + 0.2e1 * t59 * t13 + 0.2e1 * t2 * t27 + 0.2e1 * t28 * t94, 0.2e1 * t1 * t27 - 0.2e1 * t62 * t11 + 0.2e1 * t58 * t13 + 0.2e1 * t28 * t93; 0, 0, 0, 0, 0, t42, -t77, 0, -t72, t73, -t72, t56, -t73, t56 * pkin(6), 0, 0, -t12, t11, 0, t7, -t6, t5, t3, -t9, t8, 0, t54 * t48 - t98 * t51, t98 * t48 + t54 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(3) * t97, 0, 0, 0, 0, 0, 0.2e1 * t17, -0.2e1 * t16, t35, t26, 0, 0, 0, -0.2e1 * t29 * t80 + 0.2e1 * t91, -0.2e1 * t29 * t79 - 0.2e1 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t48 - t60 * t79, t55 * t51 + t60 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t81, 0, 0, 0, 0, 0, t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, -t7, t6, -t5, -t3, t9, -t8, 0, -t93 + t65 * t48 + (t13 * t48 - t64 * t51) * qJD(5), t94 + t65 * t51 + (t13 * t51 + t64 * t48) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, -0.2e1 * t68, t102, 0, 0, 0, t48 * t67 - t91, t51 * t67 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t81, 0, 0, 0, 0, 0, -t20, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t26, 0, 0, 0, t48 * t75, t51 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t59, t11, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t80, 0, t48 * t16 - t30 * t79, t51 * t16 + t30 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t81 - t49 * t79, t49 * t80 - t51 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t80, 0, -pkin(8) * t79, pkin(8) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t21;
