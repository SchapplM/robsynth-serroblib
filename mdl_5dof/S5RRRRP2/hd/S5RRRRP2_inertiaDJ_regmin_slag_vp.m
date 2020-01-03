% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:45
% EndTime: 2020-01-03 12:11:47
% DurationCPUTime: 0.45s
% Computational Cost: add. (833->99), mult. (1896->154), div. (0->0), fcn. (1570->6), ass. (0->85)
t112 = qJD(3) + qJD(4);
t73 = sin(qJ(3));
t111 = -pkin(7) - pkin(8);
t91 = qJD(3) * t111;
t51 = t73 * t91;
t76 = cos(qJ(3));
t52 = t76 * t91;
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t57 = t111 * t73;
t71 = t76 * pkin(8);
t58 = t76 * pkin(7) + t71;
t81 = t75 * t57 - t72 * t58;
t17 = -t81 * qJD(4) - t75 * t51 - t72 * t52;
t100 = pkin(1) * qJD(2);
t77 = cos(qJ(2));
t94 = t77 * t100;
t86 = t76 * t94;
t74 = sin(qJ(2));
t64 = t74 * pkin(1) + pkin(7);
t108 = -pkin(8) - t64;
t88 = qJD(3) * t108;
t37 = t73 * t88 + t86;
t87 = t73 * t94;
t38 = t76 * t88 - t87;
t46 = t108 * t73;
t47 = t76 * t64 + t71;
t84 = t75 * t46 - t72 * t47;
t10 = -t84 * qJD(4) - t75 * t37 - t72 * t38;
t110 = t75 * pkin(3);
t109 = t77 * pkin(1);
t107 = t72 * t73;
t50 = t72 * t76 + t75 * t73;
t33 = t112 * t50;
t49 = -t75 * t76 + t107;
t106 = -t33 * qJ(5) - t49 * qJD(5);
t97 = t73 * qJD(3);
t68 = pkin(3) * t97;
t69 = t74 * t100;
t53 = t69 + t68;
t67 = -t76 * pkin(3) - pkin(2);
t56 = t67 - t109;
t105 = t56 * t33 + t53 * t49;
t70 = t76 * qJD(3);
t98 = qJD(4) * t75;
t32 = t112 * t107 - t75 * t70 - t76 * t98;
t104 = -t56 * t32 + t53 * t50;
t103 = t67 * t33 + t49 * t68;
t102 = -t67 * t32 + t50 * t68;
t66 = -pkin(2) - t109;
t101 = t66 * t70 + t73 * t69;
t99 = t50 * qJ(5);
t96 = pkin(2) * t97;
t95 = pkin(2) * t70;
t93 = qJD(4) * t72 * pkin(3);
t92 = pkin(3) * t98;
t27 = t33 * pkin(4) + t68;
t19 = t84 - t99;
t45 = t49 * qJ(5);
t83 = -t72 * t46 - t75 * t47;
t20 = -t45 - t83;
t3 = -t10 + t106;
t11 = t83 * qJD(4) - t72 * t37 + t75 * t38;
t79 = t32 * qJ(5) - t50 * qJD(5);
t4 = t11 + t79;
t90 = t19 * t32 - t20 * t33 - t3 * t49 - t4 * t50;
t25 = t81 - t99;
t80 = -t72 * t57 - t75 * t58;
t26 = -t45 - t80;
t8 = t106 - t17;
t18 = t80 * qJD(4) - t72 * t51 + t75 * t52;
t9 = t18 + t79;
t89 = t25 * t32 - t26 * t33 - t8 * t49 - t9 * t50;
t78 = t66 * t97 - t76 * t69;
t40 = t49 * pkin(4) + t67;
t65 = pkin(4) + t110;
t60 = 0.2e1 * t73 * t70;
t48 = 0.2e1 * (-t73 ^ 2 + t76 ^ 2) * qJD(3);
t39 = t40 - t109;
t31 = pkin(4) * t32;
t22 = t27 + t69;
t21 = -0.2e1 * t50 * t32;
t12 = 0.2e1 * t32 * t49 - 0.2e1 * t50 * t33;
t7 = t65 * t32 + (-t33 * t72 + (-t49 * t75 + t50 * t72) * qJD(4)) * pkin(3);
t1 = [0, 0, 0, 0, -0.2e1 * t69, -0.2e1 * t94, t60, t48, 0, 0, 0, 0.2e1 * t78, 0.2e1 * t101, t21, t12, 0, 0, 0, 0.2e1 * t105, 0.2e1 * t104, 0.2e1 * t90, 0.2e1 * t19 * t4 + 0.2e1 * t20 * t3 + 0.2e1 * t39 * t22; 0, 0, 0, 0, -t69, -t94, t60, t48, 0, 0, 0, t78 - t96, -t95 + t101, t21, t12, 0, 0, 0, t103 + t105, t102 + t104, t89 + t90, t19 * t9 + t20 * t8 + t22 * t40 + t4 * t25 + t3 * t26 + t39 * t27; 0, 0, 0, 0, 0, 0, t60, t48, 0, 0, 0, -0.2e1 * t96, -0.2e1 * t95, t21, t12, 0, 0, 0, 0.2e1 * t103, 0.2e1 * t102, 0.2e1 * t89, 0.2e1 * t25 * t9 + 0.2e1 * t26 * t8 + 0.2e1 * t40 * t27; 0, 0, 0, 0, 0, 0, 0, 0, t70, -t97, 0, -t64 * t70 - t87, t64 * t97 - t86, 0, 0, -t32, -t33, 0, t11, t10, t7, t4 * t65 + (t3 * t72 + (-t19 * t72 + t20 * t75) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, t70, -t97, 0, -pkin(7) * t70, pkin(7) * t97, 0, 0, -t32, -t33, 0, t18, t17, t7, t9 * t65 + (t72 * t8 + (-t25 * t72 + t26 * t75) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t93, -0.2e1 * t92, 0, 0.2e1 * (-t65 + t110) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t11, t10, t31, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t18, t17, t31, t9 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t92, 0, -pkin(4) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
