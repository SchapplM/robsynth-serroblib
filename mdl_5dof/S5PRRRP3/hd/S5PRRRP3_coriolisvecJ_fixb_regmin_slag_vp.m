% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:15
% DurationCPUTime: 0.57s
% Computational Cost: add. (832->124), mult. (2104->179), div. (0->0), fcn. (1390->4), ass. (0->92)
t112 = cos(qJ(4));
t114 = pkin(6) + pkin(7);
t66 = cos(qJ(3));
t50 = t114 * t66;
t65 = sin(qJ(3));
t93 = t65 * qJD(1);
t34 = qJD(2) * t50 + t93;
t64 = sin(qJ(4));
t28 = t64 * t34;
t87 = qJD(2) * t114;
t33 = t66 * qJD(1) - t65 * t87;
t99 = qJD(3) * pkin(3);
t31 = t33 + t99;
t83 = t112 * t31 - t28;
t45 = t112 * t65 + t64 * t66;
t96 = qJD(2) * t45;
t97 = t96 * qJ(5);
t117 = t97 - t83;
t61 = qJD(3) + qJD(4);
t92 = qJD(2) * qJD(3);
t116 = -0.2e1 * t92;
t115 = t96 ^ 2;
t5 = t61 * pkin(4) - t117;
t113 = t5 + t117;
t107 = t64 * t65;
t77 = t61 * t107;
t81 = t112 * qJD(4);
t89 = t112 * t66;
t20 = -qJD(3) * t89 - t66 * t81 + t77;
t111 = t20 * t61;
t78 = qJD(2) * t89;
t95 = qJD(2) * t65;
t90 = t64 * t95;
t36 = -t78 + t90;
t59 = -t66 * pkin(3) - pkin(2);
t48 = t59 * qJD(2);
t22 = t36 * pkin(4) + qJD(5) + t48;
t110 = t22 * t96;
t109 = t96 * t36;
t108 = t48 * t96;
t68 = qJD(2) ^ 2;
t106 = t66 * t68;
t67 = qJD(3) ^ 2;
t105 = t67 * t65;
t104 = t67 * t66;
t21 = t61 * t45;
t18 = t21 * qJD(2);
t103 = -t45 * t18 + t20 * t36;
t102 = t112 * t33 - t28;
t101 = t61 * t78;
t100 = t65 ^ 2 - t66 ^ 2;
t98 = t36 * qJ(5);
t94 = qJD(4) * t64;
t91 = t65 * t99;
t30 = t112 * t34;
t85 = t65 * t92;
t56 = pkin(3) * t85;
t88 = t18 * pkin(4) + t56;
t86 = qJD(3) * t114;
t26 = t33 * qJD(3);
t27 = (-t66 * t87 - t93) * qJD(3);
t84 = t112 * t27 - t64 * t26;
t82 = -t64 * t33 - t30;
t80 = pkin(2) * t116;
t17 = qJD(2) * t77 - t101;
t44 = -t89 + t107;
t76 = -t44 * t17 + t21 * t96;
t75 = -t64 * t31 - t30;
t49 = t114 * t65;
t74 = -t112 * t50 + t64 * t49;
t46 = t65 * t86;
t47 = t66 * t86;
t73 = -t112 * t46 - t64 * t47 - t49 * t81 - t50 * t94;
t72 = t75 * qJD(4) + t84;
t71 = t74 * qJD(4) - t112 * t47 + t64 * t46;
t70 = t112 * t26 + t64 * t27 + t31 * t81 - t34 * t94;
t69 = t48 * t36 - t70;
t58 = t112 * pkin(3) + pkin(4);
t35 = t36 ^ 2;
t19 = t21 * t61;
t16 = -t44 * qJ(5) - t74;
t15 = -t45 * qJ(5) - t112 * t49 - t64 * t50;
t13 = -t35 + t115;
t10 = t101 + (t36 - t90) * t61;
t9 = -t97 + t102;
t8 = t82 + t98;
t7 = -t75 - t98;
t4 = t20 * qJ(5) - t45 * qJD(5) + t71;
t3 = -t21 * qJ(5) - t44 * qJD(5) + t73;
t2 = t17 * qJ(5) - qJD(5) * t96 + t72;
t1 = -t18 * qJ(5) - t36 * qJD(5) + t70;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t104, 0, 0, 0, 0, 0, -t19, t111, t76 + t103, t1 * t45 - t2 * t44 - t7 * t20 - t5 * t21; 0, 0, 0, 0, 0.2e1 * t66 * t85, t100 * t116, t104, -t105, 0, -pkin(6) * t104 + t65 * t80, pkin(6) * t105 + t66 * t80, -t17 * t45 - t20 * t96, -t76 + t103, -t111, -t19, 0, t59 * t18 + t48 * t21 + t36 * t91 + t56 * t44 + t71 * t61, -t59 * t17 - t48 * t20 - t73 * t61 + 0.2e1 * t96 * t91, -t1 * t44 + t15 * t17 - t16 * t18 - t2 * t45 + t5 * t20 - t7 * t21 - t3 * t36 - t4 * t96, t1 * t16 + t7 * t3 + t2 * t15 + t5 * t4 + t88 * (t44 * pkin(4) + t59) + t22 * (t21 * pkin(4) + t91); 0, 0, 0, 0, -t65 * t106, t100 * t68, 0, 0, 0, t68 * pkin(2) * t65, pkin(2) * t106, t109, t13, t10, 0, 0, -pkin(3) * t36 * t95 - t108 - t82 * t61 + (-t30 + (-pkin(3) * t61 - t31) * t64) * qJD(4) + t84, t102 * t61 + (-t61 * t81 - t95 * t96) * pkin(3) + t69, t58 * t17 + (t7 + t8) * t96 + (-t5 + t9) * t36 + (-t18 * t64 + (-t112 * t36 + t64 * t96) * qJD(4)) * pkin(3), -pkin(4) * t110 + t2 * t58 - t5 * t8 - t7 * t9 + (-t22 * t95 + t1 * t64 + (t112 * t7 - t5 * t64) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t13, t10, 0, 0, -t75 * t61 - t108 + t72, t83 * t61 + t69, pkin(4) * t17 - t113 * t36, t113 * t7 + (t2 - t110) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 - t115, t7 * t36 + t5 * t96 + t88;];
tauc_reg = t6;
