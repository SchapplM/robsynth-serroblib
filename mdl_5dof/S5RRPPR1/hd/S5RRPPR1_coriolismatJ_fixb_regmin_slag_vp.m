% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:56:01
% EndTime: 2020-01-03 11:56:03
% DurationCPUTime: 0.56s
% Computational Cost: add. (617->104), mult. (1383->139), div. (0->0), fcn. (1265->8), ass. (0->95)
t68 = cos(pkin(8));
t70 = sin(qJ(2));
t110 = t68 * t70;
t71 = cos(qJ(2));
t113 = t71 * pkin(1);
t61 = pkin(2) + t113;
t66 = sin(pkin(8));
t79 = pkin(1) * t110 + t66 * t61;
t37 = qJ(4) + t79;
t59 = t66 * pkin(2) + qJ(4);
t65 = sin(pkin(9));
t63 = t65 ^ 2;
t67 = cos(pkin(9));
t64 = t67 ^ 2;
t126 = (t64 / 0.2e1 + t63 / 0.2e1) * (t37 + t59);
t69 = sin(qJ(5));
t108 = t69 * t67;
t112 = cos(qJ(5));
t84 = t112 * t65;
t49 = t84 + t108;
t114 = t67 * pkin(4);
t58 = t66 * t70 * pkin(1);
t82 = t68 * t61 - t58;
t78 = -pkin(3) - t82;
t28 = t78 - t114;
t85 = -t68 * pkin(2) - pkin(3);
t52 = t85 - t114;
t87 = t52 / 0.2e1 + t28 / 0.2e1;
t124 = t87 * t49;
t109 = t69 * t65;
t83 = t112 * t67;
t47 = -t83 + t109;
t14 = t47 ^ 2 - t49 ^ 2;
t93 = qJD(1) + qJD(2);
t123 = t93 * t14;
t122 = t93 * t47;
t80 = t93 * t49;
t57 = t63 + t64;
t121 = t93 * t57;
t111 = t66 * t71;
t46 = t68 * t113 - t58;
t16 = t57 * t46;
t53 = t57 * qJD(4);
t107 = t16 * qJD(2) + t53;
t106 = pkin(1) * qJD(1);
t105 = pkin(1) * qJD(2);
t45 = (t110 + t111) * pkin(1);
t1 = t37 * t16 + t78 * t45;
t104 = t1 * qJD(1);
t103 = qJD(1) * t45;
t102 = qJD(2) * t45;
t10 = -t82 * t45 + t79 * t46;
t101 = t10 * qJD(1);
t12 = t57 * t37;
t100 = t12 * qJD(1);
t99 = t16 * qJD(1);
t98 = t47 * qJD(1);
t97 = t47 * qJD(2);
t41 = t47 * qJD(5);
t96 = t49 * qJD(1);
t95 = t49 * qJD(2);
t94 = t49 * qJD(5);
t92 = t28 * t98;
t91 = t28 * t96;
t90 = t45 * t98;
t89 = t45 * t96;
t88 = t65 * t103;
t33 = t57 * t59;
t81 = pkin(1) * t93;
t74 = (t111 / 0.2e1 + t110 / 0.2e1) * pkin(1);
t8 = t74 - t126;
t77 = t8 * qJD(1) - t33 * qJD(2);
t73 = (-t108 / 0.2e1 - t84 / 0.2e1) * t46;
t2 = t73 - t124;
t76 = -t2 * qJD(1) + t52 * t95;
t72 = (-t83 / 0.2e1 + t109 / 0.2e1) * t46;
t3 = t87 * t47 + t72;
t75 = -t3 * qJD(1) - t52 * t97;
t62 = t67 * pkin(7);
t44 = t49 * qJD(4);
t40 = t47 * qJD(4);
t39 = t67 * t59 + t62;
t38 = (-pkin(7) - t59) * t65;
t36 = t65 * t102;
t23 = t67 * t37 + t62;
t22 = (-pkin(7) - t37) * t65;
t21 = t47 * t94;
t20 = t45 * t95;
t19 = t45 * t97;
t13 = t14 * qJD(5);
t11 = t47 * t80;
t9 = t74 + t126;
t5 = t73 + t124;
t4 = t72 - (t28 + t52) * t47 / 0.2e1;
t6 = [0, 0, 0, 0, -t70 * t105, -t71 * t105, t10 * qJD(2), -t67 * t102, t36, t107, t1 * qJD(2) + t12 * qJD(4), -t21, t13, 0, 0, 0, t28 * t94 + t19, -t28 * t41 + t20; 0, 0, 0, 0, -t70 * t81, -t71 * t81, t101 + (-t45 * t68 + t46 * t66) * qJD(2) * pkin(2), -t93 * t67 * t45, t36 + t88, t99 + t107, t104 + (t46 * t33 + t45 * t85) * qJD(2) + t9 * qJD(4), -t21, t13, 0, 0, 0, t5 * qJD(5) + t19 + t90, t4 * qJD(5) + t20 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t9 * qJD(2) + t100, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t123, -t41, -t94, 0, t91 + t5 * qJD(2) + (-t112 * t23 - t69 * t22) * qJD(5), -t92 + t4 * qJD(2) + (-t112 * t22 + t69 * t23) * qJD(5); 0, 0, 0, 0, t70 * t106, t71 * t106, -t101, t67 * t103, -t88, t53 - t99, -t8 * qJD(4) - t104, -t21, t13, 0, 0, 0, -t2 * qJD(5) - t90, -t3 * qJD(5) - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t33 * qJD(4), -t21, t13, 0, 0, 0, t52 * t94, -t52 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t77, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t123, -t41, -t94, 0, (-t112 * t39 - t69 * t38) * qJD(5) + t76, (-t112 * t38 + t69 * t39) * qJD(5) + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t8 * qJD(2) - t100, 0, 0, 0, 0, 0, t94, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t77, 0, 0, 0, 0, 0, t94, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t123, 0, 0, 0, t2 * qJD(2) - t44 - t91, t3 * qJD(2) + t40 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t123, 0, 0, 0, -t44 - t76, t40 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
