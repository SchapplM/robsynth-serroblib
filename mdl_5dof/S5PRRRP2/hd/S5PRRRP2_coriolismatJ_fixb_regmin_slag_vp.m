% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP2
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
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:42:00
% EndTime: 2019-12-05 16:42:04
% DurationCPUTime: 0.66s
% Computational Cost: add. (499->121), mult. (1025->143), div. (0->0), fcn. (701->4), ass. (0->104)
t65 = cos(qJ(3));
t116 = t65 * pkin(2);
t62 = sin(qJ(4));
t64 = cos(qJ(4));
t75 = -t64 * pkin(4) - t62 * qJ(5);
t37 = -pkin(3) + t75;
t32 = t37 - t116;
t25 = t32 * t62;
t34 = t37 * t62;
t108 = t25 / 0.2e1 + t34 / 0.2e1;
t103 = t64 * qJ(5);
t117 = t62 * pkin(4);
t39 = t103 - t117;
t110 = t39 * t64;
t121 = -t108 - t110;
t54 = -pkin(3) - t116;
t86 = pkin(3) / 0.2e1 - t54 / 0.2e1;
t120 = t86 * t62;
t60 = t62 ^ 2;
t61 = t64 ^ 2;
t50 = t61 - t60;
t93 = qJD(2) + qJD(3);
t119 = t93 * t50;
t118 = pkin(3) * t64;
t115 = t32 * t39;
t114 = t32 * t64;
t113 = t37 * t64;
t112 = t39 * t37;
t111 = t39 * t62;
t109 = t54 * t64;
t104 = pkin(2) * qJD(3);
t63 = sin(qJ(3));
t92 = t63 * t104;
t49 = t62 * t92;
t57 = t60 * qJD(5);
t106 = t57 - t49;
t105 = pkin(2) * qJD(2);
t53 = t63 * pkin(2) + pkin(7);
t78 = (t60 + t61) * t65;
t7 = (t32 * t63 + t53 * t78) * pkin(2);
t102 = t7 * qJD(2);
t101 = qJD(2) * t62;
t100 = qJD(3) * t62;
t99 = qJD(4) * t62;
t58 = qJD(4) * t64;
t10 = -t111 + t114;
t98 = t10 * qJD(2);
t11 = -t25 - t110;
t97 = t11 * qJD(2);
t33 = pkin(2) * t78;
t96 = t33 * qJD(2);
t95 = t64 * qJD(5);
t94 = qJD(4) * qJ(5);
t91 = pkin(7) * t99;
t90 = pkin(7) * t58;
t89 = t63 * t105;
t88 = -t116 / 0.2e1;
t87 = t116 / 0.2e1;
t85 = t53 * t99;
t84 = t53 * t58;
t83 = qJD(2) * t115;
t82 = t32 * t101;
t81 = t54 * t101;
t80 = qJD(2) * t109;
t79 = t37 / 0.2e1 + t32 / 0.2e1;
t77 = pkin(2) * t93;
t76 = t64 * t92;
t12 = -t111 + t113;
t45 = t64 * t87;
t4 = t79 * t64 - t111 + t45;
t74 = t4 * qJD(2) + t12 * qJD(3);
t13 = -t34 - t110;
t44 = t62 * t88;
t3 = t44 + t121;
t73 = t3 * qJD(2) + t13 * qJD(3);
t51 = t62 * t95;
t72 = t51 - t76;
t71 = t39 * qJD(4) + t62 * qJD(5);
t15 = t44 + t120;
t70 = pkin(3) * t100 + t15 * qJD(2);
t46 = t64 * t88;
t16 = t86 * t64 + t46;
t69 = t16 * qJD(2) + qJD(3) * t118;
t66 = (t103 / 0.2e1 - t117 / 0.2e1) * t116;
t1 = t79 * t39 + t66;
t68 = t1 * qJD(2) + qJD(3) * t112;
t43 = t62 * t87;
t8 = t43 + t108;
t67 = t8 * qJD(2) + t37 * t100;
t26 = t75 * qJD(4) + t95;
t52 = t62 * t58;
t48 = t64 * t89;
t47 = t62 * t89;
t41 = t50 * qJD(4);
t38 = t93 * t60;
t28 = t93 * t64 * t62;
t27 = t33 * qJD(3);
t18 = -t118 / 0.2e1 + t109 / 0.2e1 + t46;
t17 = t44 - t120;
t9 = t43 - t108;
t6 = t44 - t121;
t5 = -t114 / 0.2e1 + t111 - t113 / 0.2e1 + t45;
t2 = -t112 / 0.2e1 - t115 / 0.2e1 + t66;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t58, -t99, 0, t58, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t92, -t65 * t104, t52, t41, 0, 0, 0, t54 * t99 - t76, t54 * t58 + t49, -t11 * qJD(4) + t72, t27, -t10 * qJD(4) + t106, t7 * qJD(3) - t71 * t32; 0, 0, 0, 0, 0, -t63 * t77, -t65 * t77, t52, t41, 0, 0, 0, t17 * qJD(4) - t48 - t76, t18 * qJD(4) + t47 + t49, t6 * qJD(4) - t48 + t72, t27 + t96, t5 * qJD(4) + t106 - t47, t102 + t2 * qJD(4) + t9 * qJD(5) + (pkin(7) * t78 + t37 * t63) * t104; 0, 0, 0, 0, 0, 0, 0, t28, t119, t58, -t99, 0, t17 * qJD(3) + t81 - t84, t18 * qJD(3) + t80 + t85, t6 * qJD(3) - t84 - t97, t26, t5 * qJD(3) - t85 - t98, t2 * qJD(3) + t26 * t53 - t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t58, t38, t9 * qJD(3) - t82 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t89, t65 * t105, t52, t41, 0, 0, 0, -t15 * qJD(4) + t48, -t16 * qJD(4) - t47, -t3 * qJD(4) + t48 + t51, -t96, -t4 * qJD(4) + t47 + t57, -t1 * qJD(4) - t8 * qJD(5) - t102; 0, 0, 0, 0, 0, 0, 0, t52, t41, 0, 0, 0, -pkin(3) * t99, -pkin(3) * t58, -t13 * qJD(4) + t51, 0, -t12 * qJD(4) + t57, -t71 * t37; 0, 0, 0, 0, 0, 0, 0, t28, t119, t58, -t99, 0, -t70 - t90, -t69 + t91, -t73 - t90, t26, -t74 - t91, t26 * pkin(7) - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t58, t38, -t67 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t28, -t119, 0, 0, 0, t15 * qJD(3) - t81, t16 * qJD(3) - t80, t3 * qJD(3) + t97, 0, t4 * qJD(3) + t98, t1 * qJD(3) + t83; 0, 0, 0, 0, 0, 0, 0, -t28, -t119, 0, 0, 0, t70, t69, t73, 0, t74, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t38, t8 * qJD(3) + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t38, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t14;
