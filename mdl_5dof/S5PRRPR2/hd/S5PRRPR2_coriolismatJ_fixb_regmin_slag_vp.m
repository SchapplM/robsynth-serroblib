% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:37
% EndTime: 2019-12-05 16:17:40
% DurationCPUTime: 0.67s
% Computational Cost: add. (482->113), mult. (1315->173), div. (0->0), fcn. (949->6), ass. (0->116)
t74 = sin(pkin(9));
t72 = t74 ^ 2;
t75 = cos(pkin(9));
t73 = t75 ^ 2;
t65 = t72 + t73;
t107 = qJD(2) + qJD(3);
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t46 = (t76 ^ 2 - t78 ^ 2) * t72;
t143 = t107 * t46;
t47 = t65 * t76;
t22 = t107 * t47;
t48 = t65 * t78;
t23 = t107 * t48;
t142 = t107 * t65;
t90 = t107 * t75;
t77 = sin(qJ(3));
t138 = t77 * pkin(2);
t70 = qJ(4) + t138;
t141 = qJ(4) + t70;
t140 = t70 / 0.2e1;
t139 = pkin(3) * t77;
t79 = cos(qJ(3));
t137 = t79 * pkin(2);
t136 = qJ(4) / 0.2e1;
t131 = t70 * t75;
t54 = -t75 * pkin(4) - t74 * pkin(7) - pkin(3);
t43 = t54 - t137;
t24 = t76 * t131 - t78 * t43;
t135 = t24 * t75;
t25 = -t78 * t131 - t76 * t43;
t134 = t25 * t75;
t119 = qJ(4) * t75;
t30 = t76 * t119 - t78 * t54;
t133 = t30 * t75;
t31 = -t78 * t119 - t76 * t54;
t132 = t31 * t75;
t130 = t72 * t76;
t129 = t72 * t78;
t128 = t76 * t77;
t127 = t76 * t79;
t126 = t77 * t78;
t125 = t78 * t79;
t19 = ((t75 * t125 + t128) * t75 + t72 * t125) * pkin(2);
t41 = t48 * qJD(4);
t124 = t19 * qJD(3) + t41;
t42 = t65 * t137;
t58 = t65 * qJD(4);
t123 = t42 * qJD(3) + t58;
t18 = (-t72 * t127 + (-t75 * t127 + t126) * t75) * pkin(2);
t40 = t47 * qJD(4);
t122 = -t18 * qJD(3) + t40;
t121 = pkin(2) * qJD(2);
t120 = pkin(2) * qJD(3);
t5 = -t70 * t130 - t135;
t118 = t5 * qJD(2);
t49 = t70 * t129;
t6 = -t49 + t134;
t117 = t6 * qJD(2);
t116 = t72 * qJ(4);
t115 = qJD(4) * t75;
t114 = qJD(5) * t76;
t113 = qJD(5) * t78;
t37 = t65 * t70;
t11 = (-t139 + (t37 - t138) * t79) * pkin(2);
t112 = t11 * qJD(2);
t111 = t18 * qJD(2);
t110 = t19 * qJD(2);
t109 = t37 * qJD(2);
t108 = t42 * qJD(2);
t106 = t76 * t129;
t105 = t77 * t120;
t104 = t77 * t121;
t103 = t138 / 0.2e1;
t102 = t75 * t114;
t101 = t75 * t113;
t100 = t76 * t115;
t99 = t78 * t115;
t97 = -t127 / 0.2e1;
t96 = t126 / 0.2e1;
t95 = -t125 / 0.2e1;
t66 = t78 * t116;
t94 = -t49 / 0.2e1 - t66 / 0.2e1;
t57 = t65 * qJ(4);
t92 = pkin(2) * t107;
t91 = t74 * t104;
t89 = qJD(5) * t106;
t88 = t77 * t92;
t87 = t76 * t90;
t86 = t78 * t90;
t1 = pkin(2) * t96 + (pkin(2) * t97 + t31 / 0.2e1 + t25 / 0.2e1) * t75 + t94;
t15 = -t66 + t132;
t85 = -t1 * qJD(2) - t15 * qJD(3);
t14 = -t76 * t116 - t133;
t2 = (-t138 / 0.2e1 + (t136 + t140) * t72) * t76 + (pkin(2) * t95 + t30 / 0.2e1 + t24 / 0.2e1) * t75;
t84 = -t2 * qJD(2) + t14 * qJD(3);
t9 = t103 + t141 * (-t73 / 0.2e1 - t72 / 0.2e1);
t83 = t9 * qJD(2) - t57 * qJD(3);
t82 = -qJD(5) + t90;
t81 = t82 * t76;
t80 = t82 * t78;
t62 = t74 * t105;
t56 = t74 * t101;
t55 = t74 * t102;
t39 = t46 * qJD(5);
t34 = t107 * t106;
t33 = t74 * t86;
t32 = t74 * t87;
t27 = t74 * t80;
t26 = t74 * t81;
t10 = t73 * t136 + t116 / 0.2e1 + t103 + t65 * t140;
t8 = t101 - t23;
t7 = t102 - t22;
t4 = -t132 / 0.2e1 - t134 / 0.2e1 + (t75 * t97 + t96) * pkin(2) - t94;
t3 = -t133 / 0.2e1 - t135 / 0.2e1 + (t75 * t95 - t128 / 0.2e1) * pkin(2) - t141 * t130 / 0.2e1;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t113, t74 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t105, -t79 * t120, -t75 * t105, t62, t123, t11 * qJD(3) + t37 * qJD(4), -t89, t39, t55, t56, 0, -t6 * qJD(5) + t122, t5 * qJD(5) + t124; 0, 0, 0, 0, 0, -t88, -t79 * t92, -t75 * t88, t62 + t91, t108 + t123, t112 + t10 * qJD(4) + (t79 * t57 - t139) * t120, -t89, t39, t55, t56, 0, t4 * qJD(5) - t111 + t122, t3 * qJD(5) + t110 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t10 * qJD(3) + t109, 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t143, t26, t27, 0, t4 * qJD(3) + t25 * qJD(5) - t117, t3 * qJD(3) + t24 * qJD(5) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t104, t79 * t121, t75 * t104, -t91, t58 - t108, -t9 * qJD(4) - t112, -t89, t39, t55, t56, 0, -t1 * qJD(5) + t111 + t40, -t2 * qJD(5) - t110 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57 * qJD(4), -t89, t39, t55, t56, 0, -t15 * qJD(5) + t40, t14 * qJD(5) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, -t83, 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t143, t26, t27, 0, t31 * qJD(5) + t85, t30 * qJD(5) + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t9 * qJD(3) - t109, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t83, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t143, -t32, -t33, 0, t1 * qJD(3) - t100 + t117, t2 * qJD(3) - t118 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t143, -t32, -t33, 0, -t85 - t100, -t84 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t12;
