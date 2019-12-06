% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:58
% EndTime: 2019-12-05 18:00:01
% DurationCPUTime: 0.55s
% Computational Cost: add. (991->137), mult. (2101->192), div. (0->0), fcn. (1287->4), ass. (0->100)
t65 = qJD(3) + qJD(4);
t73 = cos(qJ(3));
t105 = qJD(1) * t73;
t74 = -pkin(1) - pkin(6);
t53 = t74 * qJD(1) + qJD(2);
t31 = -pkin(7) * t105 + t73 * t53;
t27 = qJD(3) * pkin(3) + t31;
t103 = qJD(3) * t73;
t88 = pkin(7) * qJD(1) - t53;
t29 = t88 * t103;
t72 = cos(qJ(4));
t126 = (qJD(4) * t27 - t29) * t72;
t95 = t72 * t105;
t71 = sin(qJ(3));
t106 = qJD(1) * t71;
t70 = sin(qJ(4));
t96 = t70 * t106;
t39 = t95 - t96;
t32 = t39 * qJ(5);
t30 = -pkin(7) * t106 + t71 * t53;
t24 = t70 * t30;
t92 = t72 * t27 - t24;
t125 = t32 - t92;
t44 = t70 * t73 + t72 * t71;
t37 = t44 * qJD(1);
t124 = t39 ^ 2;
t66 = qJD(1) * qJD(2);
t123 = 0.2e1 * t66;
t5 = t65 * pkin(4) - t125;
t122 = t5 + t125;
t121 = pkin(7) - t74;
t102 = qJD(4) * t70;
t104 = qJD(3) * t71;
t115 = t72 * t73;
t85 = t65 * t115;
t19 = -t71 * t102 - t70 * t104 + t85;
t120 = t19 * t65;
t50 = pkin(3) * t106 + qJD(1) * qJ(2);
t20 = t37 * pkin(4) + qJD(5) + t50;
t119 = t20 * t39;
t118 = t39 * t37;
t117 = t50 * t39;
t25 = t72 * t30;
t49 = t121 * t73;
t116 = t72 * t49;
t75 = qJD(3) ^ 2;
t114 = t75 * t71;
t113 = t75 * t73;
t112 = t72 * t31 - t24;
t111 = t65 * t96;
t99 = qJD(1) * qJD(3);
t93 = t73 * t99;
t47 = pkin(3) * t93 + t66;
t110 = t71 ^ 2 - t73 ^ 2;
t76 = qJD(1) ^ 2;
t109 = -t75 - t76;
t108 = t37 * qJ(5);
t107 = t76 * qJ(2);
t59 = t71 * pkin(3) + qJ(2);
t101 = t20 * qJD(1);
t54 = pkin(3) * t103 + qJD(2);
t100 = qJ(2) * qJD(3);
t98 = 0.2e1 * qJD(1);
t97 = pkin(3) * t105;
t94 = -pkin(3) * t65 - t27;
t28 = t88 * t104;
t91 = t72 * t28 + t70 * t29;
t90 = -t30 * t102 + t70 * t28;
t89 = -t70 * t31 - t25;
t16 = qJD(1) * t85 - t111;
t86 = t16 * pkin(4) + t47;
t18 = t65 * t44;
t15 = t18 * qJD(1);
t45 = -t70 * t71 + t115;
t84 = -t45 * t15 - t39 * t18;
t83 = -t70 * t27 - t25;
t48 = t121 * t71;
t82 = t72 * t48 + t70 * t49;
t81 = t50 * t37 - t90;
t42 = t121 * t104;
t43 = qJD(3) * t49;
t80 = -qJD(4) * t116 + t48 * t102 + t70 * t42 - t72 * t43;
t1 = -t16 * qJ(5) - t37 * qJD(5) + t126 + t90;
t78 = t83 * qJD(4) + t91;
t2 = t15 * qJ(5) - t39 * qJD(5) + t78;
t7 = -t83 - t108;
t79 = t1 * t44 - t5 * t18 + t7 * t19 + t2 * t45;
t77 = t82 * qJD(4) + t72 * t42 + t70 * t43;
t60 = t72 * pkin(3) + pkin(4);
t36 = t37 ^ 2;
t17 = t18 * t65;
t14 = -t44 * qJ(5) - t82;
t13 = -t45 * qJ(5) + t70 * t48 - t116;
t12 = -t36 + t124;
t11 = t111 + (t39 - t95) * t65;
t9 = -t32 + t112;
t8 = t89 + t108;
t4 = t18 * qJ(5) - t45 * qJD(5) + t77;
t3 = -t19 * qJ(5) - t44 * qJD(5) + t80;
t6 = [0, 0, 0, 0, t123, qJ(2) * t123, -0.2e1 * t71 * t93, 0.2e1 * t110 * t99, -t114, -t113, 0, -t74 * t114 + (qJD(2) * t71 + t73 * t100) * t98, -t74 * t113 + (qJD(2) * t73 - t71 * t100) * t98, t84, t15 * t44 - t45 * t16 + t18 * t37 - t39 * t19, -t17, -t120, 0, t59 * t16 + t50 * t19 + t54 * t37 + t47 * t44 + t77 * t65, -t59 * t15 - t50 * t18 + t54 * t39 + t47 * t45 - t80 * t65, t13 * t15 - t14 * t16 - t3 * t37 - t4 * t39 - t79, t1 * t14 + t7 * t3 + t2 * t13 + t5 * t4 + t86 * (t44 * pkin(4) + t59) + t20 * (t19 * pkin(4) + t54); 0, 0, 0, 0, -t76, -t107, 0, 0, 0, 0, 0, t109 * t71, t109 * t73, 0, 0, 0, 0, 0, -qJD(1) * t37 - t17, -qJD(1) * t39 - t120, -t44 * t16 - t19 * t37 - t84, t79 - t101; 0, 0, 0, 0, 0, 0, t73 * t76 * t71, -t110 * t76, 0, 0, 0, -t73 * t107, t71 * t107, t118, t12, 0, t11, 0, -t37 * t97 - t117 - t89 * t65 + (t70 * t94 - t25) * qJD(4) + t91, -t39 * t97 + t112 * t65 + (qJD(4) * t94 + t29) * t72 + t81, t60 * t15 + (t7 + t8) * t39 + (-t5 + t9) * t37 + (-t16 * t70 + (-t37 * t72 + t39 * t70) * qJD(4)) * pkin(3), -pkin(4) * t119 + t2 * t60 - t5 * t8 - t7 * t9 + (-t73 * t101 + t1 * t70 + (-t5 * t70 + t7 * t72) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t12, 0, t11, 0, -t65 * t83 - t117 + t78, t65 * t92 - t126 + t81, pkin(4) * t15 - t122 * t37, t122 * t7 + (t2 - t119) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 - t124, t7 * t37 + t5 * t39 + t86;];
tauc_reg = t6;
