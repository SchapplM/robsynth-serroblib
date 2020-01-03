% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:30
% EndTime: 2019-12-31 19:24:33
% DurationCPUTime: 0.82s
% Computational Cost: add. (1136->174), mult. (3292->335), div. (0->0), fcn. (2800->6), ass. (0->90)
t65 = sin(qJ(2));
t66 = cos(qJ(2));
t62 = sin(pkin(5));
t90 = qJ(3) * t62;
t47 = -t66 * pkin(2) - t65 * t90 - pkin(1);
t64 = cos(pkin(5));
t75 = qJ(3) * t64 + pkin(7);
t48 = t75 * t65;
t63 = cos(pkin(8));
t104 = (t47 * t62 - t48 * t64) * t63;
t85 = qJD(3) * t65;
t27 = -t62 * t85 + (pkin(2) * t65 - t66 * t90) * qJD(2);
t70 = qJD(2) * t75;
t84 = qJD(3) * t66;
t28 = t64 * t84 - t65 * t70;
t29 = -t64 * t85 - t66 * t70;
t61 = sin(pkin(8));
t98 = t61 * t64;
t99 = t61 * t62;
t10 = t27 * t99 + t63 * t28 + t29 * t98;
t89 = qJD(2) * t65;
t5 = -(qJ(4) * t89 - qJD(4) * t66) * t62 - t10;
t103 = t27 * t63;
t94 = t64 * t66;
t39 = t65 * t61 - t63 * t94;
t86 = qJD(3) * t62;
t46 = t64 * qJD(4) + t63 * t86;
t102 = t46 * t39;
t101 = t46 * t64;
t100 = t46 * t66;
t97 = t62 * t63;
t96 = t62 * t66;
t95 = t63 * t64;
t93 = t65 * t63;
t92 = pkin(3) + qJ(5);
t49 = t75 * t66;
t43 = t61 * t49;
t91 = pkin(3) * t96 + t43;
t42 = pkin(2) * t98 + t63 * t90;
t88 = qJD(2) * t66;
t60 = t62 ^ 2;
t87 = qJD(3) * t60;
t83 = qJD(4) * t60;
t82 = qJD(4) * t61;
t40 = t61 * t94 + t93;
t81 = t40 * qJD(4);
t80 = -0.2e1 * pkin(1) * qJD(2);
t16 = t47 * t99 - t48 * t98 + t63 * t49;
t79 = t60 * t84;
t78 = t64 * t86;
t77 = -pkin(2) * t63 - pkin(3);
t76 = -qJ(4) * t61 - pkin(2);
t17 = t64 * t27 - t62 * t29;
t18 = t64 * t47 + t62 * t48;
t74 = -0.2e1 * t78;
t73 = t61 * t79;
t24 = t61 * t28;
t72 = -t29 * t95 + t24;
t30 = -t64 * qJ(4) - t42;
t69 = -t40 * qJ(4) + t18;
t13 = qJ(4) * t96 - t16;
t35 = t63 * t88 - t89 * t98;
t67 = -t35 * qJ(4) + t17 - t81;
t59 = t61 ^ 2;
t56 = t62 * t89;
t54 = t61 * t90;
t53 = t61 * t86;
t52 = t59 * t87;
t45 = -t64 * qJD(5) + t53;
t41 = pkin(2) * t95 - t54;
t38 = (-qJD(5) * t63 - t82) * t62;
t34 = (t61 * t66 + t64 * t93) * qJD(2);
t32 = (-pkin(3) * t63 + t76) * t62;
t31 = t77 * t64 + t54;
t23 = (-t92 * t63 + t76) * t62;
t22 = pkin(4) * t97 - t30;
t19 = pkin(4) * t99 + t54 + (-qJ(5) + t77) * t64;
t15 = -t43 + t104;
t14 = t91 - t104;
t12 = t39 * pkin(3) + t69;
t11 = -t39 * pkin(4) - t13;
t9 = -t24 + (t27 * t62 + t29 * t64) * t63;
t8 = t92 * t39 + t69;
t7 = t48 * t95 + t40 * pkin(4) + (qJ(5) * t66 - t47 * t63) * t62 + t91;
t6 = (-pkin(3) * t89 - t103) * t62 + t72;
t4 = t34 * pkin(3) + t67;
t3 = -t34 * pkin(4) - t5;
t2 = t35 * pkin(4) + (qJD(5) * t66 - t92 * t89 - t103) * t62 + t72;
t1 = t39 * qJD(5) + t92 * t34 + t67;
t20 = [0, 0, 0, 0.2e1 * t65 * t88, 0.2e1 * (-t65 ^ 2 + t66 ^ 2) * qJD(2), 0, 0, 0, t65 * t80, t66 * t80, 0.2e1 * t17 * t39 + 0.2e1 * t18 * t34 + 0.2e1 * (t15 * t89 - t66 * t9) * t62, 0.2e1 * t17 * t40 + 0.2e1 * t18 * t35 + 0.2e1 * (t10 * t66 - t16 * t89) * t62, -0.2e1 * t10 * t39 - 0.2e1 * t15 * t35 - 0.2e1 * t16 * t34 - 0.2e1 * t9 * t40, 0.2e1 * t16 * t10 + 0.2e1 * t15 * t9 + 0.2e1 * t18 * t17, 0.2e1 * t13 * t34 + 0.2e1 * t14 * t35 + 0.2e1 * t5 * t39 + 0.2e1 * t6 * t40, -0.2e1 * t12 * t34 - 0.2e1 * t4 * t39 + 0.2e1 * (t14 * t89 - t6 * t66) * t62, -0.2e1 * t12 * t35 - 0.2e1 * t4 * t40 + 0.2e1 * (-t13 * t89 + t5 * t66) * t62, 0.2e1 * t12 * t4 + 0.2e1 * t13 * t5 + 0.2e1 * t14 * t6, -0.2e1 * t11 * t34 + 0.2e1 * t2 * t40 - 0.2e1 * t3 * t39 + 0.2e1 * t7 * t35, -0.2e1 * t1 * t40 - 0.2e1 * t8 * t35 + 0.2e1 * (t11 * t89 - t3 * t66) * t62, 0.2e1 * t1 * t39 + 0.2e1 * t8 * t34 + 0.2e1 * (t2 * t66 - t7 * t89) * t62, 0.2e1 * t8 * t1 + 0.2e1 * t11 * t3 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, t88, -t89, 0, -pkin(7) * t88, pkin(7) * t89, t73 + t9 * t64 + (-pkin(2) * t34 - t17 * t63 + t41 * t89) * t62, t63 * t79 - t10 * t64 + (-pkin(2) * t35 + t17 * t61 - t42 * t89) * t62, -t42 * t34 - t41 * t35 + (t10 * t63 - t61 * t9 + (-t39 * t63 + t40 * t61) * qJD(3)) * t62, t10 * t42 + t9 * t41 + (-pkin(2) * t17 + (-t15 * t61 + t16 * t63) * qJD(3)) * t62, t30 * t34 + t31 * t35 - t102 + (-t5 * t63 + (qJD(3) * t40 + t6) * t61) * t62, -t73 - t32 * t34 + t6 * t64 + (t31 * t89 + t39 * t82 + t4 * t63) * t62, -t32 * t35 - t5 * t64 + (-t30 * t89 - t100 + (-t4 + t81) * t61) * t62, -t13 * t46 + t5 * t30 + t6 * t31 + t4 * t32 + (qJD(3) * t14 - qJD(4) * t12) * t99, t19 * t35 - t22 * t34 - t102 + t45 * t40 + (t2 * t61 + t3 * t63) * t62, -t23 * t35 + t3 * t64 - t38 * t40 + (-t1 * t61 + t22 * t89 - t100) * t62, -t2 * t64 + t23 * t34 + t38 * t39 + (-t1 * t63 - t19 * t89 + t45 * t66) * t62, t1 * t23 + t11 * t46 + t2 * t19 + t3 * t22 + t8 * t38 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t74, t63 * t74, 0.2e1 * t63 ^ 2 * t87 + 0.2e1 * t52, 0.2e1 * (-t41 * t61 + t42 * t63) * t86, 0.2e1 * t46 * t97 + 0.2e1 * t52, 0.2e1 * (-t63 * t83 + t78) * t61, 0.2e1 * t59 * t83 + 0.2e1 * t101, -0.2e1 * t30 * t46 + 0.2e1 * (qJD(3) * t31 - qJD(4) * t32) * t99, 0.2e1 * (t45 * t61 + t46 * t63) * t62, -0.2e1 * t38 * t99 + 0.2e1 * t101, -0.2e1 * t38 * t97 - 0.2e1 * t45 * t64, 0.2e1 * t19 * t45 + 0.2e1 * t22 * t46 + 0.2e1 * t23 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, t17, 0, -t34, -t35, t4, 0, -t35, t34, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t82, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t56, 0, t6, t35, 0, -t56, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t56, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t20;
