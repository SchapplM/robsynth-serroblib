% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:24
% DurationCPUTime: 0.83s
% Computational Cost: add. (951->153), mult. (2234->224), div. (0->0), fcn. (1610->6), ass. (0->89)
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t39 = t66 * t60 + t64 * t61;
t71 = qJD(1) * t39;
t113 = qJD(5) + t71;
t97 = t66 * t61;
t85 = qJD(1) * t97;
t94 = qJD(1) * t60;
t86 = t64 * t94;
t35 = t85 - t86;
t63 = sin(qJ(5));
t65 = cos(qJ(5));
t88 = t65 * qJD(4);
t22 = t63 * t35 - t88;
t116 = t113 * t22;
t115 = qJD(5) - t113;
t78 = t65 * t113;
t43 = qJD(4) * t86;
t28 = qJD(4) * t85 - t43;
t98 = t63 * t28;
t114 = -t113 * t78 - t98;
t111 = t71 * qJD(4);
t62 = -pkin(1) - qJ(3);
t110 = t62 * qJD(1);
t96 = t60 ^ 2 + t61 ^ 2;
t109 = t96 * qJD(3);
t44 = qJD(2) + t110;
t83 = -pkin(6) * qJD(1) + t44;
t29 = t83 * t60;
t30 = t83 * t61;
t15 = t66 * t29 + t64 * t30;
t38 = t64 * t60 - t97;
t70 = t38 * qJD(3);
t4 = -qJD(1) * t70 + t15 * qJD(4);
t107 = (t35 * pkin(4) + t113 * pkin(7)) * t113 + t4;
t75 = t64 * t29 - t66 * t30;
t10 = -qJD(4) * pkin(4) + t75;
t59 = qJD(1) * qJ(2);
t53 = qJD(3) + t59;
t42 = pkin(3) * t94 + t53;
t12 = pkin(4) * t71 - t35 * pkin(7) + t42;
t49 = t60 * pkin(3) + qJ(2);
t18 = t39 * pkin(4) + t38 * pkin(7) + t49;
t103 = -pkin(6) + t62;
t40 = t103 * t60;
t41 = t103 * t61;
t21 = t66 * t40 + t64 * t41;
t72 = t39 * qJD(3);
t3 = -qJD(1) * t72 - t75 * qJD(4);
t92 = qJD(4) * t66;
t93 = qJD(4) * t64;
t36 = -t60 * t92 - t61 * t93;
t20 = t64 * t40 - t66 * t41;
t7 = -t20 * qJD(4) - t72;
t106 = -(qJD(5) * t18 + t7) * t113 - t21 * t28 - (qJD(5) * t12 + t3) * t39 - t4 * t38 + t10 * t36;
t58 = qJD(1) * qJD(2);
t84 = 0.2e1 * t58;
t90 = qJD(5) * t63;
t5 = qJD(5) * t88 - t111 * t65 - t35 * t90;
t105 = t38 * t5;
t104 = t5 * t63;
t102 = t18 * t28;
t24 = t63 * qJD(4) + t65 * t35;
t101 = t24 * t35;
t100 = t35 * t22;
t99 = t63 * t111;
t26 = t65 * t28;
t95 = qJD(1) * t38;
t91 = qJD(5) * t38;
t37 = -t60 * t93 + t61 * t92;
t89 = t37 * qJD(4);
t79 = qJD(1) * t96;
t77 = qJD(5) * t39 + qJD(1);
t11 = qJD(4) * pkin(7) + t15;
t2 = t65 * t11 + t63 * t12;
t76 = t63 * t11 - t65 * t12;
t74 = t26 + (-t63 * t71 - t90) * t113;
t73 = t65 * t36 + t38 * t90;
t68 = -pkin(7) * t28 + (t10 - t75) * t113;
t67 = qJD(1) ^ 2;
t31 = t36 * qJD(4);
t16 = t37 * pkin(4) - t36 * pkin(7) + qJD(2);
t13 = t28 * pkin(4) + pkin(7) * t111 + t58;
t9 = t65 * t13;
t8 = t21 * qJD(4) - t70;
t6 = qJD(5) * t24 - t99;
t1 = [0, 0, 0, 0, t84, qJ(2) * t84, t60 * t84, t61 * t84, 0.2e1 * qJD(3) * t79, (t53 + t59) * qJD(2) + (-t44 - t110) * t109, t111 * t38 + t35 * t36, t111 * t39 + t38 * t28 - t35 * t37 - t36 * t71, t31, -t89, 0, 0.2e1 * t71 * qJD(2) - t8 * qJD(4) + t49 * t28 + t42 * t37, -t7 * qJD(4) - t49 * t111 + t42 * t36 + (t35 - t95) * qJD(2), -t65 * t105 + t73 * t24, (-t22 * t65 - t24 * t63) * t36 + (t104 + t6 * t65 + (-t22 * t63 + t24 * t65) * qJD(5)) * t38, t113 * t73 + t24 * t37 - t38 * t26 + t5 * t39, t38 * t98 - t22 * t37 - t6 * t39 + (-t63 * t36 + t65 * t91) * t113, t113 * t37 + t28 * t39, -t76 * t37 + t20 * t6 + t8 * t22 + t9 * t39 + (t16 * t113 + t102 + (-t10 * t38 - t11 * t39 - t113 * t21) * qJD(5)) * t65 + t106 * t63, -t2 * t37 + t20 * t5 + t8 * t24 + (-(-qJD(5) * t21 + t16) * t113 - t102 - (-qJD(5) * t11 + t13) * t39 + t10 * t91) * t63 + t106 * t65; 0, 0, 0, 0, -t67, -t67 * qJ(2), -t67 * t60, -t67 * t61, 0, (-t53 - t109) * qJD(1), 0, 0, 0, 0, 0, -qJD(1) * t71 + t31, -qJD(1) * t35 - t89, 0, 0, 0, 0, 0, -t39 * t98 - t36 * t22 + t38 * t6 + (-t37 * t63 - t77 * t65) * t113, -t39 * t26 - t36 * t24 + t105 + (-t37 * t65 + t77 * t63) * t113; 0, 0, 0, 0, 0, 0, 0, 0, -t96 * t67, t44 * t79 + t58, 0, 0, 0, 0, 0, -t43 + (t35 + t85) * qJD(4), -0.2e1 * t111, 0, 0, 0, 0, 0, t74 - t100, -t101 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t71, t35 ^ 2 - t71 ^ 2, 0, t43 + (t35 - t85) * qJD(4), 0, qJD(3) * t95 - t42 * t35, (qJD(3) + t42) * t71, t24 * t78 + t104, (t5 - t116) * t65 + (-t113 * t24 - t6) * t63, -t101 - t114, t74 + t100, -t113 * t35, -pkin(4) * t6 - t107 * t65 - t15 * t22 + t35 * t76 + t68 * t63, -pkin(4) * t5 + t107 * t63 - t15 * t24 + t2 * t35 + t68 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t22, -t22 ^ 2 + t24 ^ 2, t5 + t116, -t115 * t24 + t99, t28, -t10 * t24 - t115 * t2 - t63 * t3 + t9, t10 * t22 + t115 * t76 - t63 * t13 - t65 * t3;];
tauc_reg = t1;
