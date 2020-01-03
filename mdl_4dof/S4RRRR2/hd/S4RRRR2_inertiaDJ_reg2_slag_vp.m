% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:20
% DurationCPUTime: 0.46s
% Computational Cost: add. (616->86), mult. (1516->151), div. (0->0), fcn. (1169->6), ass. (0->76)
t103 = qJD(3) + qJD(4);
t102 = pkin(7) + pkin(6);
t64 = cos(qJ(2));
t101 = t64 * pkin(1);
t62 = sin(qJ(2));
t53 = t62 * pkin(1) + pkin(6);
t100 = -pkin(7) - t53;
t99 = cos(qJ(4));
t60 = sin(qJ(4));
t61 = sin(qJ(3));
t98 = t60 * t61;
t63 = cos(qJ(3));
t81 = t99 * t61;
t41 = t60 * t63 + t81;
t24 = t103 * t41;
t80 = t99 * t63;
t40 = -t80 + t98;
t90 = t61 * qJD(3);
t85 = pkin(3) * t90;
t92 = pkin(1) * qJD(2);
t87 = t62 * t92;
t42 = t85 + t87;
t55 = -t63 * pkin(3) - pkin(2);
t45 = t55 - t101;
t97 = t45 * t24 + t42 * t40;
t75 = t99 * qJD(4);
t23 = -qJD(3) * t80 + t103 * t98 - t63 * t75;
t96 = -t45 * t23 + t42 * t41;
t95 = t55 * t24 + t40 * t85;
t94 = -t55 * t23 + t41 * t85;
t54 = -pkin(2) - t101;
t56 = t63 * qJD(3);
t93 = t54 * t56 + t61 * t87;
t91 = qJD(4) * t60;
t89 = pkin(2) * t90;
t88 = pkin(2) * t56;
t86 = t64 * t92;
t84 = pkin(3) * t91;
t83 = t60 * t102;
t82 = t61 * t56;
t57 = t63 * pkin(7);
t37 = t63 * t53 + t57;
t67 = t100 * t81;
t21 = -t60 * t37 + t67;
t22 = t100 * t98 + t99 * t37;
t73 = t63 * t86;
t76 = qJD(3) * t100;
t65 = t61 * t76 + t73;
t74 = t61 * t86;
t66 = t63 * t76 - t74;
t4 = -qJD(4) * t67 + t37 * t91 - t60 * t66 - t99 * t65;
t5 = -t22 * qJD(4) - t60 * t65 + t99 * t66;
t79 = t21 * t23 - t22 * t24 + t4 * t40 - t5 * t41;
t58 = t61 ^ 2;
t59 = t63 ^ 2;
t78 = (t58 + t59) * t64;
t46 = t63 * pkin(6) + t57;
t70 = t102 * t99;
t69 = t61 * t70;
t11 = t103 * t69 + t46 * t91 + t83 * t56;
t72 = t61 * t83;
t28 = t99 * t46 - t72;
t12 = -t28 * qJD(4) + (-t63 * t70 + t72) * qJD(3);
t27 = -t60 * t46 - t69;
t77 = t11 * t40 - t12 * t41 + t27 * t23 - t28 * t24;
t71 = pkin(3) * t75;
t68 = t54 * t90 - t63 * t87;
t49 = -0.2e1 * t82;
t48 = 0.2e1 * t82;
t39 = 0.2e1 * (-t58 + t59) * qJD(3);
t34 = t78 * t92;
t16 = -0.2e1 * t41 * t23;
t15 = 0.2e1 * t40 * t24;
t6 = 0.2e1 * t23 * t40 - 0.2e1 * t41 * t24;
t3 = (t99 * t23 - t24 * t60 + (-t99 * t40 + t41 * t60) * qJD(4)) * pkin(3);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t87, -0.2e1 * t86, 0, 0, t48, t39, 0, t49, 0, 0, 0.2e1 * t68, 0.2e1 * t93, 0.2e1 * t34, 0.2e1 * (t53 * t78 + t54 * t62) * t92, t16, t6, 0, t15, 0, 0, 0.2e1 * t97, 0.2e1 * t96, 0.2e1 * t79, 0.2e1 * t21 * t5 - 0.2e1 * t22 * t4 + 0.2e1 * t45 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86, 0, 0, t48, t39, 0, t49, 0, 0, t68 - t89, -t88 + t93, t34, (-pkin(2) * t62 + pkin(6) * t78) * t92, t16, t6, 0, t15, 0, 0, t95 + t97, t94 + t96, t77 + t79, -t22 * t11 + t21 * t12 + t5 * t27 - t4 * t28 + t42 * t55 + t45 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t39, 0, t49, 0, 0, -0.2e1 * t89, -0.2e1 * t88, 0, 0, t16, t6, 0, t15, 0, 0, 0.2e1 * t95, 0.2e1 * t94, 0.2e1 * t77, -0.2e1 * t28 * t11 + 0.2e1 * t27 * t12 + 0.2e1 * t55 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t90, 0, -t53 * t56 - t74, t53 * t90 - t73, 0, 0, 0, 0, -t23, 0, -t24, 0, t5, t4, t3, (t99 * t5 - t4 * t60 + (-t21 * t60 + t99 * t22) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t90, 0, -pkin(6) * t56, pkin(6) * t90, 0, 0, 0, 0, -t23, 0, -t24, 0, t12, t11, t3, (t99 * t12 - t11 * t60 + (-t27 * t60 + t99 * t28) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t84, -0.2e1 * t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, 0, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, 0, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
