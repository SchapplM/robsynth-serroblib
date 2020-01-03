% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:55
% EndTime: 2019-12-31 17:59:57
% DurationCPUTime: 0.67s
% Computational Cost: add. (353->116), mult. (834->175), div. (0->0), fcn. (642->6), ass. (0->103)
t58 = sin(qJ(5));
t121 = 0.2e1 * t58;
t53 = t58 ^ 2;
t60 = cos(qJ(5));
t55 = t60 ^ 2;
t34 = t55 - t53;
t61 = cos(qJ(4));
t111 = t60 * t61;
t77 = t111 * t121;
t63 = qJD(1) * t77 - t34 * qJD(4);
t59 = sin(qJ(4));
t54 = t59 ^ 2;
t120 = t54 / 0.2e1;
t119 = t59 * pkin(7);
t118 = t61 * pkin(4);
t47 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t117 = t47 * t58;
t56 = t61 ^ 2;
t116 = t56 * t58;
t50 = t56 * t60;
t28 = t118 + t119;
t115 = t58 * t28;
t114 = t60 * t28;
t113 = t60 * t47;
t112 = t60 * t54;
t92 = t61 * qJD(4);
t43 = t58 * t92;
t102 = qJD(5) * t60;
t46 = t59 * t102;
t110 = t46 + t43;
t33 = t54 - t56;
t48 = sin(pkin(8)) * pkin(1) + qJ(3);
t74 = t59 * pkin(4) - t61 * pkin(7);
t62 = t48 + t74;
t7 = t59 * t117 - t60 * t62;
t89 = t61 * t117;
t1 = -t7 * t61 + (t89 + t114) * t59;
t109 = t1 * qJD(1);
t90 = t59 * t113;
t8 = t58 * t62 + t90;
t88 = t47 * t111;
t2 = t8 * t61 + (-t88 + t115) * t59;
t108 = t2 * qJD(1);
t3 = -t47 * t116 - t7 * t59;
t107 = t3 * qJD(1);
t4 = -t56 * t113 - t8 * t59;
t106 = t4 * qJD(1);
t105 = qJD(3) * t59;
t104 = qJD(4) * t60;
t103 = qJD(5) * t58;
t87 = 0.1e1 / 0.2e1 + t120;
t12 = (-t56 / 0.2e1 - t87) * t58;
t101 = t12 * qJD(1);
t13 = t50 / 0.2e1 + t87 * t60;
t100 = t13 * qJD(1);
t26 = t33 * t58;
t99 = t26 * qJD(1);
t27 = -t50 + t112;
t98 = t27 * qJD(1);
t97 = t33 * qJD(1);
t96 = t48 * qJD(1);
t95 = t59 * qJD(1);
t94 = t59 * qJD(4);
t93 = t61 * qJD(1);
t91 = t61 * qJD(5);
t86 = t58 * t91;
t85 = t60 * t91;
t84 = t48 * t95;
t83 = t48 * t93;
t82 = t58 * t102;
t81 = t58 * t104;
t80 = t60 * t92;
t79 = t59 * t92;
t78 = t59 * t93;
t75 = qJD(4) * t77;
t45 = t59 * t103;
t73 = t45 - t80;
t72 = (-qJD(5) - t95) * t61;
t71 = t119 / 0.2e1 + t118 / 0.2e1;
t66 = t28 / 0.2e1 + t71;
t9 = t66 * t58;
t70 = pkin(4) * t104 - t9 * qJD(1);
t10 = t66 * t60;
t69 = pkin(4) * t58 * qJD(4) + t10 * qJD(1);
t68 = t60 * t72;
t20 = (t53 / 0.2e1 - t55 / 0.2e1) * t61;
t67 = -t20 * qJD(1) + t81;
t65 = t58 * qJD(1) * t50 + t20 * qJD(4);
t25 = t34 * t56;
t64 = t25 * qJD(1) + t75;
t49 = t92 / 0.2e1;
t44 = t60 * t95;
t42 = t58 * t94;
t41 = t58 * t95;
t23 = (t95 + qJD(5) / 0.2e1) * t61;
t19 = t42 - t85;
t18 = t60 * t94 + t86;
t17 = t20 * qJD(5);
t15 = -t112 / 0.2e1 - t50 / 0.2e1 + t60 / 0.2e1;
t14 = t116 / 0.2e1 + (t120 - 0.1e1 / 0.2e1) * t58;
t6 = -t89 + t114 / 0.2e1 - t71 * t60;
t5 = -t88 - t115 / 0.2e1 + t71 * t58;
t11 = [0, 0, 0, 0, 0, qJD(3), t48 * qJD(3), -t79, t33 * qJD(4), 0, 0, 0, t48 * t92 + t105, qJD(3) * t61 - t48 * t94, -t55 * t79 - t56 * t82, -t25 * qJD(5) + t59 * t75, -t27 * qJD(4) - t59 * t86, t26 * qJD(4) - t59 * t85, t79, t1 * qJD(4) + t4 * qJD(5) + t60 * t105, -t2 * qJD(4) - t3 * qJD(5) - t58 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t96, 0, 0, 0, 0, 0, t95, t93, 0, 0, 0, 0, 0, t15 * qJD(5) + t44, t14 * qJD(5) - t41; 0, 0, 0, 0, 0, 0, 0, -t78, t97, -t94, -t92, 0, -t47 * t94 + t83, -t47 * t92 - t84, -t17 + (-t55 * t93 - t81) * t59, t63 * t59 - 0.2e1 * t61 * t82, t43 - t98, t80 + t99, t23, t109 + (t58 * t74 - t90) * qJD(4) + t6 * qJD(5), -t108 + (-pkin(7) * t111 + (pkin(4) * t60 + t117) * t59) * qJD(4) + t5 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t64, t58 * t72, t68, t49, t15 * qJD(3) + t6 * qJD(4) - t8 * qJD(5) + t106, t14 * qJD(3) + t5 * qJD(4) + t7 * qJD(5) - t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t94, 0, 0, 0, 0, 0, t73, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, -qJD(1), -t96, 0, 0, 0, 0, 0, -t95, -t93, 0, 0, 0, 0, 0, -t13 * qJD(5) - t44, -t12 * qJD(5) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t92, 0, 0, 0, 0, 0, -t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 - t110, t73 - t101; 0, 0, 0, 0, 0, 0, 0, t78, -t97, 0, 0, 0, -t83, t84, t55 * t78 - t17, t68 * t121, t46 + t98, -t45 - t99, -t23, -t10 * qJD(5) - t109, t9 * qJD(5) + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t34 * qJD(5), 0, 0, 0, -pkin(4) * t103, -pkin(4) * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t63, t44 + t102, -t41 - t103, -t93 / 0.2e1, -pkin(7) * t102 - t69, pkin(7) * t103 - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, (t58 * t93 - t104) * t59, t60 * t78 + t42, t49, t13 * qJD(3) + t10 * qJD(4) - t106, t12 * qJD(3) - t9 * qJD(4) + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t63, -t44, t41, t93 / 0.2e1, t69, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
