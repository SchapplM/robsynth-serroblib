% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:39
% EndTime: 2019-03-09 02:18:42
% DurationCPUTime: 0.83s
% Computational Cost: add. (1610->118), mult. (3502->201), div. (0->0), fcn. (3657->10), ass. (0->87)
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t57 = sin(pkin(11));
t58 = cos(pkin(11));
t112 = t105 * t58 + t107 * t57;
t113 = t112 * qJD(3);
t70 = t105 * t57 - t107 * t58;
t68 = t70 * qJD(4);
t61 = cos(qJ(6));
t56 = t61 ^ 2;
t60 = sin(qJ(6));
t96 = t60 ^ 2 - t56;
t82 = t96 * qJD(6);
t111 = qJD(4) + qJD(5);
t49 = sin(pkin(10)) * pkin(1) + qJ(3);
t109 = pkin(7) + t49;
t38 = t109 * t57;
t39 = t109 * t58;
t71 = t105 * t39 + t107 * t38;
t110 = t70 * qJD(3) + t71 * qJD(4);
t104 = sin(qJ(5));
t106 = cos(qJ(5));
t27 = -pkin(8) * t112 - t71;
t72 = t105 * t38 - t107 * t39;
t28 = -t70 * pkin(8) - t72;
t16 = t104 * t27 + t106 * t28;
t85 = qJD(4) * t105;
t86 = qJD(4) * t107;
t62 = -pkin(8) * t68 - t38 * t85 + t39 * t86 + t113;
t97 = t57 * t86 + t58 * t85;
t63 = -t97 * pkin(8) - t110;
t5 = t16 * qJD(5) + t104 * t63 + t106 * t62;
t3 = t5 * t60;
t15 = t104 * t28 - t106 * t27;
t52 = qJD(6) * t61;
t108 = t15 * t52 + t3;
t66 = t106 * t70;
t83 = qJD(5) * t104;
t24 = t104 * t97 + t111 * t66 + t112 * t83;
t65 = t104 * t70;
t31 = t106 * t112 - t65;
t103 = t31 * t24;
t84 = qJD(5) * t106;
t25 = t106 * t97 - t111 * t65 + t112 * t84;
t102 = t60 * t25;
t101 = t61 * t24;
t100 = t61 * t25;
t30 = t104 * t112 + t66;
t99 = t31 * t100 - t30 * t101;
t51 = -t106 * pkin(4) - pkin(5);
t80 = pkin(4) * t83;
t98 = t51 * t52 + t60 * t80;
t95 = qJD(6) * t60;
t94 = t60 * t101;
t93 = pkin(5) * t95;
t92 = pkin(5) * t52;
t91 = t31 * t95;
t90 = t60 * t52;
t89 = t97 * pkin(4);
t81 = pkin(4) * t84;
t79 = 0.2e1 * (t57 ^ 2 + t58 ^ 2) * qJD(3);
t42 = -cos(pkin(10)) * pkin(1) - t58 * pkin(3) - pkin(2);
t32 = t70 * pkin(4) + t42;
t19 = t30 * pkin(5) - t31 * pkin(9) + t32;
t78 = t61 * t16 + t60 * t19;
t77 = t60 * t16 - t61 * t19;
t76 = t24 * t30 - t31 * t25;
t50 = t104 * pkin(4) + pkin(9);
t75 = t30 * t50 - t31 * t51;
t74 = -t60 * t24 + t31 * t52;
t10 = t91 + t101;
t11 = t30 * t52 + t102;
t9 = t30 * t95 - t100;
t73 = t51 * t95 - t61 * t80;
t67 = -0.2e1 * t68;
t64 = -t24 * t51 - t25 * t50 + (t104 * t31 - t106 * t30) * qJD(5) * pkin(4);
t46 = 0.2e1 * t90;
t41 = -0.2e1 * t82;
t29 = t31 ^ 2;
t13 = t15 * t95;
t8 = -t31 * t82 - t94;
t7 = t25 * pkin(5) + t24 * pkin(9) + t89;
t6 = t96 * t24 - 0.4e1 * t31 * t90;
t4 = t104 * t62 - t106 * t63 - t27 * t84 + t28 * t83;
t2 = -t78 * qJD(6) + t60 * t4 + t61 * t7;
t1 = t77 * qJD(6) + t61 * t4 - t60 * t7;
t12 = [0, 0, 0, 0, 0, 0, t79, t49 * t79, t112 * t67, 0.2e1 * t70 ^ 2 * qJD(4) - 0.2e1 * t112 * t97, 0, 0, 0, 0.2e1 * t42 * t97, t42 * t67, -0.2e1 * t103, 0.2e1 * t76, 0, 0, 0, 0.2e1 * t32 * t25 + 0.2e1 * t30 * t89, -0.2e1 * t32 * t24 + 0.2e1 * t31 * t89, -0.2e1 * t56 * t103 - 0.2e1 * t29 * t90, 0.2e1 * t29 * t82 + 0.4e1 * t31 * t94, -0.2e1 * t30 * t91 + 0.2e1 * t99, -0.2e1 * t31 * t102 - 0.2e1 * t74 * t30, 0.2e1 * t30 * t25, 0.2e1 * t74 * t15 + 0.2e1 * t2 * t30 - 0.2e1 * t77 * t25 + 0.2e1 * t31 * t3, 0.2e1 * t5 * t61 * t31 + 0.2e1 * t1 * t30 - 0.2e1 * t10 * t15 - 0.2e1 * t78 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 * t61 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -t68, 0, 0, 0, 0, 0, t25, -t24, 0, 0, 0, 0, 0, -t9, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t97, 0, t72 * qJD(4) - t113, t110, 0, 0, -t24, -t25, 0, -t5, t4, t8, t6, t11, -t9, 0, t13 + (-t75 * qJD(6) - t5) * t61 + t64 * t60, t64 * t61 + t75 * t95 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t68, 0, 0, 0, 0, 0, -t25, t24, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t80, -0.2e1 * t81, t46, t41, 0, 0, 0, 0.2e1 * t73, 0.2e1 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, -t5, t4, t8, t6, t11, -t9, 0, t13 + (pkin(5) * t24 - pkin(9) * t25) * t60 + (-t5 + (-pkin(5) * t31 - pkin(9) * t30) * qJD(6)) * t61, t10 * pkin(5) + t9 * pkin(9) + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t81, t46, t41, 0, 0, 0, t73 - t93, -t92 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t41, 0, 0, 0, -0.2e1 * t93, -0.2e1 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t74, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t95, 0, -t50 * t52 - t60 * t81, t50 * t95 - t61 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t95, 0, -pkin(9) * t52, pkin(9) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
