% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:02
% DurationCPUTime: 0.69s
% Computational Cost: add. (437->108), mult. (1177->184), div. (0->0), fcn. (870->8), ass. (0->85)
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t50 = sin(pkin(9));
t96 = qJD(1) * t50;
t110 = t57 * qJD(2) - t54 * t96;
t56 = cos(qJ(4));
t46 = -t56 * pkin(4) - pkin(3);
t18 = t46 * qJD(3) - t110;
t114 = t18 + t110;
t105 = pkin(6) + pkin(7);
t29 = t54 * qJD(2) + t57 * t96;
t47 = qJD(4) + qJD(5);
t106 = qJD(5) - t47;
t52 = sin(qJ(5));
t53 = sin(qJ(4));
t55 = cos(qJ(5));
t31 = t52 * t56 + t55 * t53;
t112 = t31 * t47;
t10 = qJD(3) * t112;
t84 = qJD(3) * qJD(4);
t113 = -0.2e1 * t84;
t30 = t52 * t53 - t55 * t56;
t62 = t30 * t47;
t58 = qJD(4) ^ 2;
t59 = qJD(3) ^ 2;
t111 = (t58 + t59) * t54;
t90 = qJD(4) * t56;
t108 = -qJD(5) * t56 - t90;
t91 = qJD(4) * t53;
t107 = qJD(5) * t53 + t91;
t93 = qJD(3) * t56;
t78 = t55 * t93;
t94 = qJD(3) * t53;
t79 = t52 * t94;
t25 = -t78 + t79;
t27 = -t52 * t93 - t55 * t94;
t104 = t27 * t25;
t103 = t50 * t57;
t51 = cos(pkin(9));
t95 = qJD(1) * t51;
t40 = t53 * t95;
t73 = t105 * qJD(3) + t29;
t12 = t73 * t56 - t40;
t102 = t55 * t12;
t101 = t59 * t54;
t100 = t59 * t57;
t99 = t53 ^ 2 - t56 ^ 2;
t97 = qJD(3) * pkin(3);
t92 = qJD(3) * t57;
t83 = t50 * t100;
t82 = pkin(4) * t94;
t81 = pkin(4) * t91;
t80 = qJD(3) * t50 * t54;
t11 = -t73 * t53 - t56 * t95;
t8 = qJD(4) * pkin(4) + t11;
t76 = -pkin(4) * t47 - t8;
t75 = qJD(4) * t105;
t74 = t56 * t84;
t19 = -t110 - t97;
t22 = t110 * qJD(3);
t72 = -t19 * qJD(3) - t22;
t71 = t53 * t80;
t70 = t56 * t80;
t69 = t57 * t113;
t68 = -t29 + t81;
t2 = t11 * qJD(4) + t56 * t22;
t3 = t40 * qJD(4) - t53 * t22 - t73 * t90;
t66 = t18 * t27 - t52 * t2 + t55 * t3;
t24 = t56 * t103 - t51 * t53;
t23 = -t53 * t103 - t51 * t56;
t65 = pkin(6) * t58;
t64 = qJD(4) * (t19 + t110 - t97);
t9 = qJD(5) * t78 - t47 * t79 + t55 * t74;
t60 = t18 * t25 + (t106 * t12 - t3) * t52;
t35 = t105 * t56;
t34 = t105 * t53;
t33 = t56 * t75;
t32 = t53 * t75;
t17 = (t29 + t81) * qJD(3);
t16 = t23 * qJD(4) - t70;
t15 = -t24 * qJD(4) + t71;
t6 = -t25 ^ 2 + t27 ^ 2;
t5 = -t27 * t47 - t10;
t4 = t25 * t47 + t9;
t1 = [0, 0, 0, -t83, t50 * t101, 0, 0, 0, 0, 0, -t56 * t83 + (t15 + t71) * qJD(4), t53 * t83 + (-t16 + t70) * qJD(4), 0, 0, 0, 0, 0, (t55 * t15 - t52 * t16 + (-t23 * t52 - t24 * t55) * qJD(5)) * t47 + (t54 * t10 + t25 * t92) * t50, -(t52 * t15 + t55 * t16 + (t23 * t55 - t24 * t52) * qJD(5)) * t47 + (-t27 * t92 + t54 * t9) * t50; 0, 0, 0, -t101, -t100, 0, 0, 0, 0, 0, -t56 * t111 + t53 * t69, t53 * t111 + t56 * t69, 0, 0, 0, 0, 0, -0.2e1 * t57 * t10 + ((t107 * t52 + t108 * t55) * t47 + qJD(3) * t25) * t54, (qJD(3) * t62 - t9) * t57 + (-(-t107 * t55 + t108 * t52) * t47 - qJD(3) * t27) * t54; 0, 0, 0, 0, 0, 0.2e1 * t53 * t74, t99 * t113, t58 * t56, -t58 * t53, 0, t53 * t64 - t65 * t56, t65 * t53 + t56 * t64, t27 * t62 + t9 * t31, -t31 * t10 + t112 * t27 + t25 * t62 - t9 * t30, -t62 * t47, -t112 * t47, 0, (t52 * t32 - t55 * t33 + (t34 * t52 - t35 * t55) * qJD(5)) * t47 + t46 * t10 + t17 * t30 + t68 * t25 + t114 * t112, -(-t55 * t32 - t52 * t33 + (-t34 * t55 - t35 * t52) * qJD(5)) * t47 + t46 * t9 + t17 * t31 - t68 * t27 - t114 * t62; 0, 0, 0, 0, 0, -t53 * t59 * t56, t99 * t59, 0, 0, 0, t72 * t53, t72 * t56, -t104, t6, t4, t5, 0, -(-t52 * t11 - t102) * t47 - t25 * t82 + (t76 * t52 - t102) * qJD(5) + t66, t27 * t82 + (t76 * qJD(5) + t11 * t47 - t2) * t55 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t6, t4, t5, 0, t66 + t106 * (-t52 * t8 - t102), (-t106 * t8 - t2) * t55 + t60;];
tauc_reg = t1;
