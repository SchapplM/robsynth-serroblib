% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRR5
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
% MMD_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:15
% DurationCPUTime: 0.59s
% Computational Cost: add. (511->115), mult. (1455->235), div. (0->0), fcn. (1183->6), ass. (0->88)
t47 = sin(qJ(2));
t101 = -0.4e1 * t47;
t45 = sin(qJ(4));
t46 = sin(qJ(3));
t48 = cos(qJ(4));
t49 = cos(qJ(3));
t25 = t45 * t49 + t48 * t46;
t17 = t25 * t47;
t50 = cos(qJ(2));
t58 = -t50 * pkin(2) - t47 * pkin(6);
t30 = -pkin(1) + t58;
t91 = t49 * t50;
t37 = pkin(5) * t91;
t88 = t46 * t30 + t37;
t81 = t50 * qJD(2);
t66 = t49 * t81;
t85 = qJD(3) * t46;
t100 = -t47 * t85 + t66;
t42 = t47 ^ 2;
t61 = (-t50 ^ 2 + t42) * qJD(2);
t43 = t49 ^ 2;
t87 = t46 ^ 2 - t43;
t62 = t87 * qJD(3);
t99 = qJD(3) + qJD(4);
t57 = pkin(2) * t47 - pkin(6) * t50;
t28 = t57 * qJD(2);
t40 = t47 * qJD(2);
t83 = qJD(3) * t50;
t72 = t46 * t83;
t52 = t49 * t40 + t72;
t84 = qJD(3) * t49;
t10 = pkin(5) * t52 - t46 * t28 - t30 * t84;
t98 = pkin(6) + pkin(7);
t97 = pkin(5) * t46;
t69 = t46 * t81;
t51 = t47 * t84 + t69;
t5 = -pkin(7) * t51 - t10;
t96 = t48 * t5;
t94 = t46 * t47;
t15 = -pkin(7) * t94 + t88;
t95 = t45 * t15;
t93 = t47 * t49;
t92 = t48 * t15;
t70 = t46 * t40;
t89 = pkin(5) * t70 + t49 * t28;
t82 = qJD(4) * t45;
t80 = -0.2e1 * pkin(1) * qJD(2);
t79 = -0.2e1 * pkin(2) * qJD(3);
t78 = pkin(3) * t85;
t77 = pkin(3) * t40;
t76 = pkin(3) * t82;
t75 = qJD(4) * t48 * pkin(3);
t74 = pkin(5) * t81;
t71 = t49 * t83;
t68 = t46 * t84;
t67 = t47 * t81;
t4 = (pkin(3) * t47 - pkin(7) * t91) * qJD(2) + (-t37 + (pkin(7) * t47 - t30) * t46) * qJD(3) + t89;
t65 = t48 * t4 - t45 * t5;
t23 = t49 * t30;
t12 = -pkin(7) * t93 + t23 + (-pkin(3) - t97) * t50;
t64 = t50 * pkin(3) - t12;
t63 = qJD(3) * t98;
t60 = 0.2e1 * t67;
t59 = t46 * t66;
t56 = t48 * t12 - t95;
t55 = t45 * t12 + t92;
t33 = t98 * t46;
t34 = t98 * t49;
t54 = -t48 * t33 - t45 * t34;
t53 = -t45 * t33 + t48 * t34;
t24 = t45 * t46 - t48 * t49;
t39 = -t49 * pkin(3) - pkin(2);
t36 = -0.2e1 * t67;
t29 = (pkin(3) * t46 + pkin(5)) * t47;
t27 = t49 * t63;
t26 = t46 * t63;
t18 = t24 * t47;
t16 = pkin(3) * t51 + t74;
t14 = t99 * t25;
t13 = t99 * t24;
t11 = -t88 * qJD(3) + t89;
t9 = -qJD(4) * t53 + t45 * t26 - t48 * t27;
t8 = -qJD(4) * t54 + t48 * t26 + t45 * t27;
t7 = -t82 * t94 + (t99 * t93 + t69) * t48 + t100 * t45;
t6 = -t99 * t17 - t24 * t81;
t2 = -qJD(4) * t55 + t65;
t1 = -qJD(4) * t56 - t45 * t4 - t96;
t3 = [0, 0, 0, t60, -0.2e1 * t61, 0, 0, 0, t47 * t80, t50 * t80, -0.2e1 * t42 * t68 + 0.2e1 * t43 * t67, t59 * t101 + 0.2e1 * t42 * t62, 0.2e1 * t47 * t72 + 0.2e1 * t49 * t61, -0.2e1 * t46 * t61 + 0.2e1 * t47 * t71, t36, 0.2e1 * t23 * t40 - 0.2e1 * t11 * t50 + 0.2e1 * (t42 * t84 + t46 * t67) * pkin(5), -0.2e1 * t10 * t50 - 0.2e1 * t88 * t40 + 0.2e1 * (-t42 * t85 + t49 * t60) * pkin(5), -0.2e1 * t18 * t6, -0.2e1 * t6 * t17 + 0.2e1 * t18 * t7, -0.2e1 * t18 * t40 - 0.2e1 * t6 * t50, -0.2e1 * t17 * t40 + 0.2e1 * t7 * t50, t36, 0.2e1 * t16 * t17 - 0.2e1 * t2 * t50 + 0.2e1 * t29 * t7 + 0.2e1 * t56 * t40, -0.2e1 * t1 * t50 - 0.2e1 * t16 * t18 + 0.2e1 * t29 * t6 - 0.2e1 * t55 * t40; 0, 0, 0, 0, 0, t81, -t40, 0, -t74, pkin(5) * t40, -t47 * t62 + t59, t68 * t101 - t87 * t81, t70 - t71, t52, 0, (pkin(6) * t91 + (-pkin(2) * t49 + t97) * t47) * qJD(3) + (t46 * t58 - t37) * qJD(2), (pkin(5) * t93 + t46 * t57) * qJD(3) + (t49 * t58 + t50 * t97) * qJD(2), t18 * t13 + t6 * t25, t13 * t17 + t18 * t14 - t6 * t24 - t25 * t7, t13 * t50 + t25 * t40, t14 * t50 - t24 * t40, 0, t29 * t14 + t16 * t24 + t17 * t78 + t39 * t7 + t54 * t40 - t9 * t50, -t29 * t13 + t16 * t25 - t18 * t78 + t39 * t6 - t53 * t40 - t8 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t68, -0.2e1 * t62, 0, 0, 0, t46 * t79, t49 * t79, -0.2e1 * t25 * t13, 0.2e1 * t13 * t24 - 0.2e1 * t25 * t14, 0, 0, 0, 0.2e1 * t39 * t14 + 0.2e1 * t24 * t78, -0.2e1 * t39 * t13 + 0.2e1 * t25 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t51, t40, t11, t10, 0, 0, t6, -t7, t40, t48 * t77 + (t45 * t64 - t92) * qJD(4) + t65, -t96 + (-t4 - t77) * t45 + (t48 * t64 + t95) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t85, 0, -pkin(6) * t84, pkin(6) * t85, 0, 0, -t13, -t14, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76, -0.2e1 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7, t40, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
