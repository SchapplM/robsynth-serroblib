% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:26
% EndTime: 2019-03-09 02:57:28
% DurationCPUTime: 0.56s
% Computational Cost: add. (770->107), mult. (1584->190), div. (0->0), fcn. (1411->6), ass. (0->76)
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t80 = sin(pkin(9));
t70 = qJD(3) * t80;
t81 = cos(pkin(9));
t71 = qJD(3) * t81;
t24 = t49 * t70 - t51 * t71;
t29 = -t81 * t49 - t80 * t51;
t89 = t29 * t24;
t25 = -t49 * t71 - t51 * t70;
t28 = -t80 * t49 + t81 * t51;
t90 = t28 * t25;
t95 = 0.2e1 * t89 + 0.2e1 * t90;
t94 = (-t80 * t24 + t81 * t25) * pkin(3);
t48 = sin(qJ(6));
t46 = t48 ^ 2;
t50 = cos(qJ(6));
t84 = -t50 ^ 2 + t46;
t72 = t84 * qJD(6);
t93 = 2 * qJD(2);
t92 = 2 * qJD(5);
t91 = -pkin(4) - pkin(8);
t88 = t29 * t50;
t87 = t48 * t25;
t86 = t50 * t24;
t85 = t50 * t25;
t83 = t49 * pkin(3) + qJ(2);
t52 = -pkin(1) - pkin(7);
t82 = qJ(4) - t52;
t79 = qJD(6) * t48;
t78 = qJD(6) * t50;
t77 = t49 * qJD(3);
t76 = t51 * qJD(3);
t34 = pkin(3) * t76 + qJD(2);
t75 = qJ(2) * qJD(3);
t74 = t48 * t78;
t73 = 0.4e1 * t48 * t88;
t33 = t82 * t51;
t27 = t29 ^ 2;
t69 = qJD(6) * (t28 ^ 2 + t27);
t68 = -t28 * qJ(5) + t83;
t43 = -t81 * pkin(3) - pkin(4);
t23 = -qJD(3) * t33 - t49 * qJD(4);
t57 = -t51 * qJD(4) + t82 * t77;
t10 = t81 * t23 + t80 * t57;
t32 = t82 * t49;
t20 = -t80 * t32 + t81 * t33;
t21 = -t81 * t32 - t80 * t33;
t9 = t80 * t23 - t81 * t57;
t67 = t21 * t10 + t20 * t9;
t11 = t28 * pkin(5) + t20;
t8 = t91 * t29 + t68;
t66 = t50 * t11 - t48 * t8;
t65 = t48 * t11 + t50 * t8;
t39 = t80 * pkin(3) + qJ(5);
t63 = t29 * qJD(5) + t24 * t39;
t17 = -t48 * t24 - t29 * t78;
t14 = t29 * t79 - t86;
t62 = -t28 * t78 - t87;
t15 = t28 * t79 - t85;
t60 = -t25 * qJ(5) - t28 * qJD(5) + t34;
t59 = t43 * t25 + t63;
t38 = -pkin(8) + t43;
t6 = t24 * pkin(5) + t10;
t58 = t6 + (-t28 * t38 - t29 * t39) * qJD(6);
t55 = t10 * t29 + t20 * t25 + t21 * t24 + t9 * t28;
t12 = t29 * pkin(5) + t21;
t54 = qJD(6) * t12 + t25 * t38 + t63;
t53 = 0.2e1 * t55;
t19 = -t29 * pkin(4) + t68;
t7 = -t24 * pkin(4) + t60;
t5 = t25 * pkin(5) + t9;
t4 = t91 * t24 + t60;
t2 = -t65 * qJD(6) - t48 * t4 + t50 * t5;
t1 = -t66 * qJD(6) - t50 * t4 - t48 * t5;
t3 = [0, 0, 0, 0, t93, qJ(2) * t93, -0.2e1 * t49 * t76, 0.2e1 * (t49 ^ 2 - t51 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t49 + 0.2e1 * t51 * t75, 0.2e1 * qJD(2) * t51 - 0.2e1 * t49 * t75, t53, 0.2e1 * t83 * t34 + 0.2e1 * t67, t53, 0.2e1 * t19 * t24 + 0.2e1 * t7 * t29, -0.2e1 * t19 * t25 - 0.2e1 * t7 * t28, 0.2e1 * t19 * t7 + 0.2e1 * t67, 0.2e1 * t27 * t74 + 0.2e1 * t46 * t89, t24 * t73 - 0.2e1 * t27 * t72, 0.2e1 * t17 * t28 - 0.2e1 * t29 * t87, 0.2e1 * t14 * t28 - 0.2e1 * t29 * t85, 0.2e1 * t90, -0.2e1 * t14 * t12 + 0.2e1 * t2 * t28 + 0.2e1 * t66 * t25 + 0.2e1 * t6 * t88, -0.2e1 * t6 * t48 * t29 + 0.2e1 * t1 * t28 + 0.2e1 * t17 * t12 - 0.2e1 * t65 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t55, -t95, 0, 0, -t55, 0, 0, 0, 0, 0, t48 * t69 - t50 * t95, t48 * t95 + t50 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, -t52 * t77, -t52 * t76, -t94 (t80 * t10 - t81 * t9) * pkin(3), t59, t9, t10, t21 * qJD(5) + t10 * t39 + t9 * t43, t29 * t72 - t48 * t86, qJD(6) * t73 + t84 * t24, -t15, t62, 0, t58 * t48 + t54 * t50, -t54 * t48 + t58 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, t94, 0, -t25, -t24, -t59, 0, 0, 0, 0, 0, t17, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t39 * t92, -0.2e1 * t74, 0.2e1 * t72, 0, 0, 0, 0.2e1 * qJD(5) * t48 + 0.2e1 * t39 * t78, 0.2e1 * qJD(5) * t50 - 0.2e1 * t39 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t24, -t25, t7, 0, 0, 0, 0, 0, t62, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, t9, 0, 0, 0, 0, 0, -t15, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t14, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t78, 0, -t38 * t79, -t38 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
