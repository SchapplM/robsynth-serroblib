% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:28
% EndTime: 2019-03-08 18:47:31
% DurationCPUTime: 1.30s
% Computational Cost: add. (1170->190), mult. (3729->375), div. (0->0), fcn. (3925->14), ass. (0->107)
t68 = sin(pkin(13));
t72 = cos(pkin(13));
t76 = sin(qJ(6));
t79 = cos(qJ(6));
t120 = -t76 * t68 + t79 * t72;
t64 = -t72 * pkin(5) - pkin(4);
t119 = 0.2e1 * t64;
t77 = sin(qJ(4));
t118 = t68 * t77;
t80 = cos(qJ(4));
t117 = t68 * t80;
t70 = sin(pkin(7));
t78 = sin(qJ(3));
t116 = t70 * t78;
t81 = cos(qJ(3));
t115 = t70 * t81;
t114 = t72 * t77;
t113 = t72 * t80;
t73 = cos(pkin(12));
t74 = cos(pkin(7));
t112 = t73 * t74;
t109 = pkin(10) + qJ(5);
t41 = -t77 * qJD(5) + (pkin(4) * t77 - qJ(5) * t80) * qJD(4);
t105 = t77 * qJD(4);
t99 = pkin(9) * t105;
t28 = t72 * t41 + t68 * t99;
t92 = -t80 * pkin(4) - t77 * qJ(5);
t55 = -pkin(3) + t92;
t62 = pkin(9) * t113;
t36 = t68 * t55 + t62;
t108 = qJD(3) * t78;
t107 = qJD(5) * t80;
t106 = qJD(6) * t77;
t104 = t80 * qJD(4);
t103 = pkin(9) * t117;
t102 = -0.2e1 * pkin(3) * qJD(4);
t101 = t81 * t112;
t100 = t77 * t116;
t65 = pkin(9) * t104;
t98 = t68 * t104;
t97 = t70 * t108;
t96 = qJD(3) * t115;
t95 = t77 * t104;
t94 = 0.2e1 * (t68 ^ 2 + t72 ^ 2) * qJD(5);
t69 = sin(pkin(12));
t71 = sin(pkin(6));
t75 = cos(pkin(6));
t27 = t75 * t116 + (t78 * t112 + t69 * t81) * t71;
t42 = -t71 * t73 * t70 + t75 * t74;
t13 = t27 * t77 - t42 * t80;
t23 = -t75 * t96 + (-qJD(3) * t101 + t108 * t69) * t71;
t10 = -t13 * qJD(4) - t23 * t80;
t24 = t27 * qJD(3);
t7 = -t10 * t68 + t24 * t72;
t8 = t10 * t72 + t24 * t68;
t93 = -t7 * t68 + t8 * t72;
t14 = t27 * t80 + t42 * t77;
t26 = -t75 * t115 + (t69 * t78 - t101) * t71;
t11 = -t14 * t68 + t26 * t72;
t12 = t14 * t72 + t26 * t68;
t91 = t79 * t11 - t76 * t12;
t90 = t76 * t11 + t79 * t12;
t32 = qJD(4) * t100 - t74 * t104 - t80 * t96;
t18 = t68 * t32 + t72 * t97;
t19 = -t72 * t32 + t68 * t97;
t89 = -t18 * t68 + t19 * t72;
t49 = t72 * t55;
t25 = -pkin(10) * t114 + t49 + (-pkin(9) * t68 - pkin(5)) * t80;
t34 = -pkin(10) * t118 + t36;
t88 = t79 * t25 - t76 * t34;
t87 = t76 * t25 + t79 * t34;
t37 = t68 * t41;
t29 = -t72 * t99 + t37;
t86 = -t28 * t68 + t29 * t72;
t46 = t80 * t116 + t77 * t74;
t30 = -t72 * t115 - t68 * t46;
t31 = -t68 * t115 + t72 * t46;
t85 = t79 * t30 - t76 * t31;
t84 = t76 * t30 + t79 * t31;
t58 = t109 * t68;
t59 = t109 * t72;
t83 = -t79 * t58 - t76 * t59;
t82 = -t76 * t58 + t79 * t59;
t52 = t79 * t68 + t76 * t72;
t53 = (pkin(5) * t68 + pkin(9)) * t77;
t47 = pkin(5) * t98 + t65;
t45 = -t80 * t74 + t100;
t44 = t52 * qJD(6);
t43 = t120 * qJD(6);
t40 = t120 * t77;
t39 = t52 * t77;
t35 = t49 - t103;
t33 = t46 * qJD(4) + t77 * t96;
t22 = t37 + (-pkin(9) * t114 - pkin(10) * t117) * qJD(4);
t21 = t52 * t104 + t120 * t106;
t20 = t104 * t120 - t52 * t106;
t17 = (pkin(5) * t77 - pkin(10) * t113) * qJD(4) + t28;
t16 = -t52 * qJD(5) - t82 * qJD(6);
t15 = -qJD(5) * t120 - t83 * qJD(6);
t9 = t14 * qJD(4) - t23 * t77;
t6 = -t87 * qJD(6) + t79 * t17 - t76 * t22;
t5 = -t88 * qJD(6) - t76 * t17 - t79 * t22;
t4 = -t84 * qJD(6) + t79 * t18 - t76 * t19;
t3 = -t85 * qJD(6) - t76 * t18 - t79 * t19;
t2 = -t90 * qJD(6) + t79 * t7 - t76 * t8;
t1 = -t91 * qJD(6) - t76 * t7 - t79 * t8;
t38 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t11 * t7 + 0.2e1 * t12 * t8 + 0.2e1 * t13 * t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t18 + t12 * t19 + t13 * t33 + t7 * t30 + t8 * t31 + t9 * t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t30 * t18 + 0.2e1 * t31 * t19 + 0.2e1 * t45 * t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t24, t23, 0, 0, 0, 0, 0, t26 * t105 - t24 * t80, t26 * t104 + t24 * t77, t9 * t118 - t7 * t80 + (t11 * t77 + t13 * t117) * qJD(4), t9 * t114 + t8 * t80 + (t13 * t113 - t12 * t77) * qJD(4) (-t68 * t8 - t7 * t72) * t77 + (-t11 * t72 - t12 * t68) * t104, t11 * t28 + t12 * t29 + t7 * t35 + t8 * t36 + (t104 * t13 + t77 * t9) * pkin(9), 0, 0, 0, 0, 0, t105 * t91 + t13 * t21 - t2 * t80 + t9 * t39, -t1 * t80 - t105 * t90 + t13 * t20 + t9 * t40; 0, 0, 0, -t97, -t96, 0, 0, 0, 0, 0 (-t81 * t105 - t80 * t108) * t70 (-t81 * t104 + t77 * t108) * t70, t33 * t118 - t18 * t80 + (t45 * t117 + t30 * t77) * qJD(4), t33 * t114 + t19 * t80 + (t45 * t113 - t31 * t77) * qJD(4) (-t18 * t72 - t19 * t68) * t77 + (-t30 * t72 - t31 * t68) * t104, t18 * t35 + t19 * t36 + t30 * t28 + t31 * t29 + (t104 * t45 + t33 * t77) * pkin(9), 0, 0, 0, 0, 0, t105 * t85 + t45 * t21 + t33 * t39 - t4 * t80, -t105 * t84 + t45 * t20 - t3 * t80 + t33 * t40; 0, 0, 0, 0, 0, 0.2e1 * t95, 0.2e1 * (-t77 ^ 2 + t80 ^ 2) * qJD(4), 0, 0, 0, t77 * t102, t80 * t102, -0.2e1 * t28 * t80 + 0.2e1 * (t35 + 0.2e1 * t103) * t105, 0.2e1 * t29 * t80 + 0.2e1 * (-t36 + 0.2e1 * t62) * t105, 0.2e1 * (-t28 * t72 - t29 * t68) * t77 + 0.2e1 * (-t35 * t72 - t36 * t68) * t104, 0.2e1 * pkin(9) ^ 2 * t95 + 0.2e1 * t35 * t28 + 0.2e1 * t36 * t29, 0.2e1 * t40 * t20, -0.2e1 * t20 * t39 - 0.2e1 * t40 * t21, 0.2e1 * t105 * t40 - 0.2e1 * t20 * t80, -0.2e1 * t105 * t39 + 0.2e1 * t21 * t80, -0.2e1 * t95, 0.2e1 * t105 * t88 + 0.2e1 * t53 * t21 + 0.2e1 * t47 * t39 - 0.2e1 * t6 * t80, -0.2e1 * t105 * t87 + 0.2e1 * t53 * t20 + 0.2e1 * t47 * t40 - 0.2e1 * t5 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -t9 * t72, t9 * t68, t93, -t9 * pkin(4) + (-t11 * t68 + t12 * t72) * qJD(5) + t93 * qJ(5), 0, 0, 0, 0, 0, -t120 * t9 + t13 * t44, t13 * t43 + t9 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, -t33 * t72, t33 * t68, t89, -t33 * pkin(4) + (-t30 * t68 + t31 * t72) * qJD(5) + t89 * qJ(5), 0, 0, 0, 0, 0, -t120 * t33 + t45 * t44, t33 * t52 + t45 * t43; 0, 0, 0, 0, 0, 0, 0, t104, -t105, 0, -t65, t99, t68 * t107 + (t68 * t92 - t62) * qJD(4), t72 * t107 + (t72 * t92 + t103) * qJD(4), t86, -pkin(4) * t65 + (-t35 * t68 + t36 * t72) * qJD(5) + t86 * qJ(5), t20 * t52 + t40 * t43, t120 * t20 - t52 * t21 - t43 * t39 - t40 * t44, t105 * t52 - t43 * t80, t105 * t120 + t44 * t80, 0, t105 * t83 - t120 * t47 - t16 * t80 + t64 * t21 + t53 * t44, -t105 * t82 - t15 * t80 + t64 * t20 + t53 * t43 + t47 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, qJ(5) * t94, 0.2e1 * t52 * t43, 0.2e1 * t120 * t43 - 0.2e1 * t52 * t44, 0, 0, 0, t44 * t119, t43 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t72 * t104, 0, t65, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t21, t105, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t38;
