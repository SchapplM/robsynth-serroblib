% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:51
% EndTime: 2019-12-31 19:27:52
% DurationCPUTime: 0.60s
% Computational Cost: add. (543->109), mult. (945->135), div. (0->0), fcn. (753->6), ass. (0->89)
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t70 = sin(qJ(2));
t72 = cos(qJ(2));
t38 = (t67 * t70 + t68 * t72) * pkin(1);
t107 = -t38 / 0.2e1;
t103 = t72 * pkin(1);
t82 = -pkin(2) - t103;
t55 = -pkin(3) + t82;
t104 = t70 * pkin(1);
t56 = qJ(3) + t104;
t19 = t68 * t55 - t67 * t56;
t16 = pkin(4) - t19;
t73 = -pkin(2) - pkin(3);
t42 = -t67 * qJ(3) + t68 * t73;
t39 = pkin(4) - t42;
t116 = t16 / 0.2e1 + t39 / 0.2e1;
t115 = t107 - t116;
t65 = qJD(1) + qJD(2);
t71 = cos(qJ(5));
t110 = t65 * t71;
t35 = t67 * t110;
t69 = sin(qJ(5));
t111 = t65 * t69;
t34 = t67 * t111;
t114 = t68 * t110;
t113 = t68 * t111;
t52 = -t69 ^ 2 + t71 ^ 2;
t112 = t65 * t52;
t61 = t67 * qJD(3);
t37 = t67 * t103 - t68 * t104;
t89 = t37 * qJD(1);
t109 = t61 - t89;
t30 = t37 * qJD(2);
t50 = t71 * t61;
t102 = t71 * t30 + t50;
t101 = t30 + t61;
t62 = t68 * qJD(3);
t100 = t38 * qJD(2) + t62;
t20 = t67 * t55 + t68 * t56;
t43 = t68 * qJ(3) + t67 * t73;
t99 = pkin(1) * qJD(1);
t98 = pkin(1) * qJD(2);
t3 = -t19 * t37 + t20 * t38;
t97 = t3 * qJD(1);
t77 = t107 + t116;
t4 = t77 * t69;
t96 = t4 * qJD(1);
t5 = t77 * t71;
t95 = t5 * qJD(1);
t8 = -t19 * t67 + t20 * t68;
t94 = t8 * qJD(1);
t93 = qJD(1) * t16;
t92 = qJD(2) * t39;
t64 = qJD(5) * t69;
t91 = qJD(5) * t71;
t18 = (t56 * t72 + t82 * t70) * pkin(1);
t90 = t18 * qJD(1);
t88 = t38 * qJD(1);
t60 = t72 * t98;
t87 = t60 + qJD(3);
t86 = t70 * t98;
t85 = t71 * t89;
t84 = -t42 / 0.2e1 - t19 / 0.2e1;
t83 = t43 / 0.2e1 + t20 / 0.2e1;
t2 = (t37 / 0.2e1 + t83) * t68 + (t107 + t84) * t67;
t9 = -t67 * t42 + t68 * t43;
t76 = t2 * qJD(1) + t9 * qJD(2);
t75 = -t62 + t93;
t74 = -t62 + t92;
t66 = qJ(3) * qJD(3);
t59 = t72 * t99;
t58 = t70 * t99;
t53 = t69 * t91;
t51 = t65 * qJ(3);
t49 = t56 * qJD(3);
t48 = t52 * qJD(5);
t45 = t65 * t68;
t44 = t65 * t67;
t41 = -t58 - t86;
t40 = -pkin(7) + t43;
t36 = t69 * t110;
t17 = -pkin(7) + t20;
t15 = t68 * t64 - t35;
t14 = t68 * t91 + t34;
t7 = t71 * t115;
t6 = t69 * t115;
t1 = (-t37 / 0.2e1 + t83) * t68 + (t38 / 0.2e1 + t84) * t67;
t10 = [0, 0, 0, 0, -t86, -t60, -t86, t87, t18 * qJD(2) + t49, t101, t100, t3 * qJD(2) + t8 * qJD(3), t53, t48, 0, 0, 0, -t16 * t64 + t102, -t101 * t69 - t16 * t91; 0, 0, 0, 0, t41, -t59 - t60, t41, t59 + t87, t90 + t49 + (-pkin(2) * t70 + qJ(3) * t72) * t98, t89 + t101, t88 + t100, t97 + (-t37 * t42 + t38 * t43) * qJD(2) + t1 * qJD(3), t53, t48, 0, 0, 0, t6 * qJD(5) + t102 + t85, t7 * qJD(5) + (-t65 * t37 - t61) * t69; 0, 0, 0, 0, 0, 0, 0, t65, t65 * t56, t44, t45, t1 * qJD(2) + t94, 0, 0, 0, 0, 0, t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t112, -t91, t64, 0, t6 * qJD(2) - t17 * t91 - t69 * t93, t7 * qJD(2) + t17 * t64 - t71 * t93; 0, 0, 0, 0, t58, t59, t58, -t59 + qJD(3), t66 - t90, t109, t62 - t88, t2 * qJD(3) - t97, t53, t48, 0, 0, 0, -t4 * qJD(5) + t50 - t85, -t5 * qJD(5) - t109 * t69; 0, 0, 0, 0, 0, 0, 0, qJD(3), t66, t61, t62, t9 * qJD(3), t53, t48, 0, 0, 0, -t39 * t64 + t50, -t39 * t91 - t69 * t61; 0, 0, 0, 0, 0, 0, 0, t65, t51, t44, t45, t76, 0, 0, 0, 0, 0, t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t112, -t91, t64, 0, -t40 * t91 - t69 * t92 - t96, t40 * t64 - t71 * t92 - t95; 0, 0, 0, 0, 0, 0, 0, -t65, -qJ(3) * qJD(2) - t56 * qJD(1), -t44, -t45, -t2 * qJD(2) - t94, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, -t65, -t51, -t44, -t45, -t76, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t91 + t113, t67 * t64 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t112, 0, 0, 0, t4 * qJD(2) + t75 * t69, t5 * qJD(2) + t75 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t112, 0, 0, 0, t74 * t69 + t96, t74 * t71 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
