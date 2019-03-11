% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:47
% EndTime: 2019-03-09 02:33:51
% DurationCPUTime: 1.08s
% Computational Cost: add. (1669->116), mult. (3473->210), div. (0->0), fcn. (3651->8), ass. (0->84)
t107 = cos(qJ(5));
t59 = sin(pkin(10));
t60 = cos(pkin(10));
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t80 = t64 * t59 - t66 * t60;
t113 = t80 * qJD(4);
t40 = t66 * t59 + t64 * t60;
t63 = sin(qJ(5));
t32 = -t107 * t80 - t63 * t40;
t39 = t40 * qJD(4);
t18 = -t32 * qJD(5) + t107 * t113 + t63 * t39;
t65 = cos(qJ(6));
t101 = t65 * t18;
t75 = t107 * t40 - t63 * t80;
t62 = sin(qJ(6));
t96 = qJD(6) * t62;
t76 = t75 * t96 + t101;
t19 = -t75 * qJD(5) - t107 * t39 + t113 * t63;
t102 = t65 * t19;
t12 = -t32 * t96 + t102;
t104 = t62 * t18;
t53 = qJD(6) * t65;
t120 = -t53 * t75 + t104;
t116 = t18 * t75;
t106 = t32 * t19;
t78 = t62 * t19 + t32 * t53;
t28 = t32 ^ 2;
t61 = -pkin(1) - qJ(3);
t109 = -pkin(7) + t61;
t43 = t109 * t59;
t44 = t109 * t60;
t81 = t64 * t43 - t66 * t44;
t26 = t80 * pkin(8) - t81;
t74 = t40 * qJD(3);
t114 = t26 * qJD(4) - t74;
t58 = t65 ^ 2;
t98 = t62 ^ 2 - t58;
t88 = t98 * qJD(6);
t45 = (t59 ^ 2 + t60 ^ 2) * qJD(3);
t110 = 2 * qJD(2);
t82 = -t66 * t43 - t64 * t44;
t27 = -t40 * pkin(8) - t82;
t16 = t107 * t27 + t63 * t26;
t68 = t80 * qJD(3) + t82 * qJD(4);
t67 = -t39 * pkin(8) - t68;
t5 = t16 * qJD(5) + t107 * t67 + t114 * t63;
t3 = t5 * t62;
t15 = -t107 * t26 + t63 * t27;
t108 = t15 * t53 + t3;
t105 = t32 * t65;
t52 = -t107 * pkin(4) - pkin(5);
t97 = qJD(5) * t63;
t92 = pkin(4) * t97;
t100 = t52 * t53 + t62 * t92;
t50 = t59 * pkin(3) + qJ(2);
t95 = qJ(2) * qJD(2);
t94 = pkin(5) * t96;
t93 = pkin(5) * t53;
t91 = t62 * t53;
t90 = -0.4e1 * t62 * t105;
t89 = qJD(5) * t107;
t87 = pkin(4) * t89;
t33 = t40 * pkin(4) + t50;
t17 = pkin(5) * t75 - t32 * pkin(9) + t33;
t86 = t65 * t16 + t62 * t17;
t85 = t62 * t16 - t65 * t17;
t84 = -t75 ^ 2 - t28;
t51 = t63 * pkin(4) + pkin(9);
t83 = -t32 * t52 + t51 * t75;
t79 = t52 * t96 - t65 * t92;
t34 = -pkin(4) * t113 + qJD(2);
t71 = -0.2e1 * t106 + 0.2e1 * t116;
t69 = t19 * t52 + t18 * t51 + (-t107 * t75 + t32 * t63) * qJD(5) * pkin(4);
t49 = 0.2e1 * t91;
t42 = -0.2e1 * t88;
t13 = t15 * t96;
t8 = t62 * t102 - t32 * t88;
t7 = -pkin(5) * t18 - t19 * pkin(9) + t34;
t6 = qJD(6) * t90 - t98 * t19;
t4 = -t107 * t114 - t26 * t89 + t27 * t97 + t63 * t67;
t2 = -t86 * qJD(6) + t62 * t4 + t65 * t7;
t1 = t85 * qJD(6) + t65 * t4 - t62 * t7;
t9 = [0, 0, 0, 0, t110, 0.2e1 * t95, t59 * t110, t60 * t110, 0.2e1 * t45, -0.2e1 * t61 * t45 + 0.2e1 * t95, 0.2e1 * t80 * t39, -0.2e1 * t113 * t80 + 0.2e1 * t39 * t40, 0, 0, 0, 0.2e1 * qJD(2) * t40 - 0.2e1 * t113 * t50, -0.2e1 * qJD(2) * t80 - 0.2e1 * t50 * t39, 0.2e1 * t106, 0.2e1 * t18 * t32 - 0.2e1 * t19 * t75, 0, 0, 0, -0.2e1 * t18 * t33 + 0.2e1 * t34 * t75, 0.2e1 * t33 * t19 + 0.2e1 * t34 * t32, 0.2e1 * t58 * t106 - 0.2e1 * t28 * t91, t19 * t90 + 0.2e1 * t28 * t88, -0.2e1 * t32 * t101 + 0.2e1 * t12 * t75, 0.2e1 * t32 * t104 - 0.2e1 * t75 * t78, -0.2e1 * t116, 0.2e1 * t78 * t15 + 0.2e1 * t18 * t85 + 0.2e1 * t2 * t75 + 0.2e1 * t32 * t3, 0.2e1 * t1 * t75 + 0.2e1 * t5 * t105 + 0.2e1 * t12 * t15 + 0.2e1 * t18 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t53 + t71 * t62, t71 * t65 - t84 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, -t113, -t39, 0, 0, 0, 0, 0, -t18, t19, 0, 0, 0, 0, 0, -t76, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t113, 0, t68, t81 * qJD(4) + t74, 0, 0, t19, t18, 0, -t5, t4, t8, t6, -t120, -t76, 0, t13 + (-t83 * qJD(6) - t5) * t65 + t69 * t62, t69 * t65 + t83 * t96 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t113, 0, 0, 0, 0, 0, t19, t18, 0, 0, 0, 0, 0, t12, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, -0.2e1 * t87, t49, t42, 0, 0, 0, 0.2e1 * t79, 0.2e1 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18, 0, -t5, t4, t8, t6, -t120, -t76, 0, t13 + (-pkin(5) * t19 + pkin(9) * t18) * t62 + (-t5 + (-pkin(5) * t32 - pkin(9) * t75) * qJD(6)) * t65, -pkin(5) * t12 + t76 * pkin(9) + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18, 0, 0, 0, 0, 0, t12, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t87, t49, t42, 0, 0, 0, t79 - t94, -t93 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t42, 0, 0, 0, -0.2e1 * t94, -0.2e1 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t78, -t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t96, 0, -t51 * t53 - t62 * t87, t51 * t96 - t65 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t96, 0, -pkin(9) * t53, pkin(9) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
