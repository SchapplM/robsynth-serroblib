% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:09
% EndTime: 2019-12-05 16:10:11
% DurationCPUTime: 0.61s
% Computational Cost: add. (918->146), mult. (2353->215), div. (0->0), fcn. (1584->6), ass. (0->93)
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t118 = -qJ(4) - pkin(6);
t93 = qJD(3) * t118;
t42 = t76 * qJD(4) + t74 * t93;
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t83 = -t74 * qJD(4) + t76 * t93;
t121 = t73 * t76;
t53 = t72 * t74 - t121;
t77 = cos(qJ(2));
t86 = t53 * t77;
t116 = qJD(1) * t86 + t73 * t42 + t72 * t83;
t103 = (qJD(2) * qJD(3));
t130 = -2 * t103;
t54 = t72 * t76 + t73 * t74;
t48 = t54 * qJD(2);
t44 = t48 ^ 2;
t112 = qJD(2) * t74;
t99 = qJD(2) * t121;
t45 = t72 * t112 - t99;
t129 = -t45 ^ 2 - t44;
t75 = sin(qJ(2));
t78 = qJD(3) ^ 2;
t79 = qJD(2) ^ 2;
t128 = (t78 + t79) * t75;
t47 = t54 * qJD(3);
t108 = qJD(3) * t76;
t109 = qJD(3) * t74;
t50 = t73 * t108 - t72 * t109;
t107 = t75 * qJD(1);
t127 = pkin(3) * t109 - t107;
t104 = qJ(4) * qJD(3);
t61 = qJD(2) * pkin(6) + t107;
t106 = t77 * qJD(1);
t89 = qJD(4) + t106;
t25 = -t61 * t109 + (-t74 * t104 + t76 * t89) * qJD(2);
t80 = -t61 * t108 + (-t76 * t104 - t74 * t89) * qJD(2);
t2 = t72 * t25 - t73 * t80;
t59 = t118 * t76;
t97 = t118 * t74;
t27 = -t72 * t59 - t73 * t97;
t126 = t2 * t27;
t38 = t54 * t75;
t125 = t2 * t38;
t91 = qJ(4) * qJD(2) + t61;
t41 = t91 * t76;
t124 = t72 * t41;
t30 = t73 * t41;
t120 = t78 * t74;
t119 = t78 * t76;
t117 = -t54 * t106 + t72 * t42 - t73 * t83;
t3 = t73 * t25 + t72 * t80;
t40 = t91 * t74;
t32 = qJD(3) * pkin(3) - t40;
t12 = t72 * t32 + t30;
t96 = t74 * t103;
t52 = pkin(3) * t96 + qJD(2) * t107;
t115 = t74 ^ 2 - t76 ^ 2;
t113 = qJD(2) * pkin(2);
t111 = qJD(2) * t75;
t15 = -t73 * t40 - t124;
t105 = qJD(5) - t15;
t98 = -t76 * pkin(3) - pkin(2);
t95 = t76 * t103;
t94 = -t47 * pkin(4) + t50 * qJ(5) + t54 * qJD(5) - t127;
t90 = t77 * t130;
t11 = t73 * t32 - t124;
t88 = qJD(2) * t113;
t51 = t98 * qJD(2) + qJD(4) - t106;
t9 = t45 * pkin(4) - t48 * qJ(5) + t51;
t87 = t9 * t48 + t2;
t16 = -t77 * t48 - t50 * t75;
t17 = -qJD(2) * t86 - t75 * t47;
t36 = qJD(2) * t47;
t60 = t72 * t96;
t37 = t73 * t95 - t60;
t39 = t53 * t75;
t85 = -t16 * t48 - t17 * t45 + t39 * t36 + t38 * t37;
t84 = t36 * pkin(4) - t37 * qJ(5) + t52;
t82 = -0.2e1 * qJD(3) * t113;
t28 = -t73 * t59 + t72 * t97;
t81 = -t116 * t45 + t117 * t48 + t2 * t54 + t27 * t37 - t28 * t36;
t68 = -t73 * pkin(3) - pkin(4);
t66 = t72 * pkin(3) + qJ(5);
t26 = t53 * pkin(4) - t54 * qJ(5) + t98;
t18 = pkin(3) * t112 + t48 * pkin(4) + t45 * qJ(5);
t14 = -t72 * t40 + t30;
t10 = qJD(3) * qJ(5) + t12;
t8 = -qJD(3) * pkin(4) + qJD(5) - t11;
t4 = -t48 * qJD(5) + t84;
t1 = qJD(3) * qJD(5) + t3;
t5 = [0, 0, -t79 * t75, -t79 * t77, 0, 0, 0, 0, 0, -t76 * t128 + t74 * t90, t74 * t128 + t76 * t90, t85, t11 * t16 + t51 * t111 + t12 * t17 - t3 * t39 - t52 * t77 + t125, t16 * qJD(3) + t45 * t111 - t77 * t36, t85, t17 * qJD(3) - t48 * t111 + t77 * t37, -t1 * t39 + t10 * t17 + t9 * t111 - t8 * t16 - t4 * t77 + t125; 0, 0, 0, 0, 0.2e1 * t74 * t95, t115 * t130, t119, -t120, 0, -pkin(6) * t119 + t74 * t82, pkin(6) * t120 + t76 * t82, -t11 * t50 - t12 * t47 - t3 * t53 + t81, -t117 * t11 + t116 * t12 + t127 * t51 + t3 * t28 + t52 * t98 + t126, -t117 * qJD(3) + t26 * t36 + t4 * t53 - t94 * t45 + t9 * t47, -t1 * t53 - t10 * t47 + t8 * t50 + t81, t116 * qJD(3) - t26 * t37 - t4 * t54 + t94 * t48 - t9 * t50, t1 * t28 + t116 * t10 + t117 * t8 + t4 * t26 - t94 * t9 + t126; 0, 0, 0, 0, -t74 * t79 * t76, t115 * t79, 0, 0, 0, t74 * t88, t76 * t88, (t12 - t14) * t48 + (-t11 + t15) * t45 + (-t36 * t72 - t37 * t73) * pkin(3), t11 * t14 - t12 * t15 + (-t51 * t112 - t2 * t73 + t3 * t72) * pkin(3), t14 * qJD(3) - t18 * t45 - t87, -t66 * t36 + t68 * t37 + (t10 - t14) * t48 + (t8 - t105) * t45, t18 * t48 - t9 * t45 + (0.2e1 * qJD(5) - t15) * qJD(3) + t3, t1 * t66 + t105 * t10 - t8 * t14 - t9 * t18 + t2 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t11 * t48 + t12 * t45 + t52, 0.2e1 * t48 * qJD(3), t129, t60 + (t45 - t99) * qJD(3), t10 * t45 + (-qJD(5) - t8) * t48 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t45, -t60 + (t45 + t99) * qJD(3), -t44 - t78, -t10 * qJD(3) + t87;];
tauc_reg = t5;
