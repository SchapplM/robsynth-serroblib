% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:16
% EndTime: 2019-12-31 17:14:18
% DurationCPUTime: 0.51s
% Computational Cost: add. (494->117), mult. (1016->143), div. (0->0), fcn. (693->4), ass. (0->104)
t64 = cos(qJ(2));
t115 = t64 * pkin(1);
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t74 = -t63 * pkin(3) - t61 * qJ(4);
t37 = -pkin(2) + t74;
t32 = t37 - t115;
t25 = t32 * t61;
t34 = t37 * t61;
t107 = t25 / 0.2e1 + t34 / 0.2e1;
t102 = t63 * qJ(4);
t116 = t61 * pkin(3);
t39 = -t102 + t116;
t109 = t39 * t63;
t120 = t109 - t107;
t54 = -pkin(2) - t115;
t85 = pkin(2) / 0.2e1 - t54 / 0.2e1;
t119 = t85 * t61;
t59 = t61 ^ 2;
t60 = t63 ^ 2;
t50 = t60 - t59;
t92 = qJD(1) + qJD(2);
t118 = t92 * t50;
t117 = pkin(2) * t63;
t114 = t32 * t39;
t113 = t32 * t63;
t112 = t37 * t63;
t111 = t39 * t37;
t110 = t39 * t61;
t108 = t54 * t63;
t103 = pkin(1) * qJD(2);
t62 = sin(qJ(2));
t91 = t62 * t103;
t49 = t61 * t91;
t57 = t59 * qJD(4);
t105 = t57 - t49;
t104 = pkin(1) * qJD(1);
t53 = t62 * pkin(1) + pkin(6);
t77 = (t59 + t60) * t64;
t7 = (t32 * t62 + t53 * t77) * pkin(1);
t101 = t7 * qJD(1);
t100 = qJD(1) * t61;
t99 = qJD(2) * t61;
t98 = qJD(3) * t61;
t58 = qJD(3) * t63;
t10 = t110 + t113;
t97 = t10 * qJD(1);
t11 = -t25 + t109;
t96 = t11 * qJD(1);
t33 = pkin(1) * t77;
t95 = t33 * qJD(1);
t94 = t63 * qJD(4);
t93 = qJD(3) * qJ(4);
t90 = pkin(6) * t98;
t89 = pkin(6) * t58;
t88 = t62 * t104;
t87 = -t115 / 0.2e1;
t86 = t115 / 0.2e1;
t84 = t53 * t98;
t83 = t53 * t58;
t82 = qJD(1) * t114;
t81 = t32 * t100;
t80 = t54 * t100;
t79 = qJD(1) * t108;
t78 = t37 / 0.2e1 + t32 / 0.2e1;
t76 = pkin(1) * t92;
t75 = t63 * t91;
t12 = t110 + t112;
t45 = t63 * t86;
t4 = t78 * t63 + t110 + t45;
t73 = t4 * qJD(1) + t12 * qJD(2);
t13 = -t34 + t109;
t44 = t61 * t87;
t3 = t44 + t120;
t72 = t3 * qJD(1) + t13 * qJD(2);
t51 = t61 * t94;
t71 = t51 - t75;
t70 = qJD(3) * t39 - qJD(4) * t61;
t15 = t44 + t119;
t69 = pkin(2) * t99 + t15 * qJD(1);
t46 = t63 * t87;
t16 = t85 * t63 + t46;
t68 = t16 * qJD(1) + qJD(2) * t117;
t65 = (t102 / 0.2e1 - t116 / 0.2e1) * t115;
t1 = -t78 * t39 + t65;
t67 = t1 * qJD(1) - qJD(2) * t111;
t43 = t61 * t86;
t8 = t43 + t107;
t66 = t8 * qJD(1) + t37 * t99;
t26 = t74 * qJD(3) + t94;
t52 = t61 * t58;
t48 = t63 * t88;
t47 = t61 * t88;
t41 = t50 * qJD(3);
t38 = t92 * t59;
t28 = t92 * t63 * t61;
t27 = t33 * qJD(2);
t18 = -t117 / 0.2e1 + t108 / 0.2e1 + t46;
t17 = t44 - t119;
t9 = t43 - t107;
t6 = t44 - t120;
t5 = -t113 / 0.2e1 - t110 - t112 / 0.2e1 + t45;
t2 = t111 / 0.2e1 + t114 / 0.2e1 + t65;
t14 = [0, 0, 0, 0, -t91, -t64 * t103, t52, t41, 0, 0, 0, t54 * t98 - t75, t54 * t58 + t49, -t11 * qJD(3) + t71, t27, -t10 * qJD(3) + t105, t7 * qJD(2) + t32 * t70; 0, 0, 0, 0, -t62 * t76, -t64 * t76, t52, t41, 0, 0, 0, t17 * qJD(3) - t48 - t75, t18 * qJD(3) + t47 + t49, t6 * qJD(3) - t48 + t71, t27 + t95, t5 * qJD(3) + t105 - t47, t101 + t2 * qJD(3) + t9 * qJD(4) + (pkin(6) * t77 + t37 * t62) * t103; 0, 0, 0, 0, 0, 0, t28, t118, t58, -t98, 0, t17 * qJD(2) + t80 - t83, t18 * qJD(2) + t79 + t84, t6 * qJD(2) - t83 - t96, t26, t5 * qJD(2) - t84 - t97, t2 * qJD(2) + t26 * t53 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t58, t38, t9 * qJD(2) - t81 + t83; 0, 0, 0, 0, t88, t64 * t104, t52, t41, 0, 0, 0, -t15 * qJD(3) + t48, -t16 * qJD(3) - t47, -t3 * qJD(3) + t48 + t51, -t95, -t4 * qJD(3) + t47 + t57, -t1 * qJD(3) - t8 * qJD(4) - t101; 0, 0, 0, 0, 0, 0, t52, t41, 0, 0, 0, -pkin(2) * t98, -pkin(2) * t58, -t13 * qJD(3) + t51, 0, -t12 * qJD(3) + t57, t70 * t37; 0, 0, 0, 0, 0, 0, t28, t118, t58, -t98, 0, -t69 - t89, -t68 + t90, -t72 - t89, t26, -t73 - t90, pkin(6) * t26 - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t58, t38, -t66 + t89; 0, 0, 0, 0, 0, 0, -t28, -t118, 0, 0, 0, t15 * qJD(2) - t80, t16 * qJD(2) - t79, t3 * qJD(2) + t96, 0, t4 * qJD(2) + t97, t1 * qJD(2) - t82; 0, 0, 0, 0, 0, 0, -t28, -t118, 0, 0, 0, t69, t68, t72, 0, t73, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t38, t8 * qJD(2) + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t38, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t14;
