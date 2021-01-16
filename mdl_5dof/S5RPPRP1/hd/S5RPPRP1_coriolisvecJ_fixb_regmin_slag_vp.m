% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:55:47
% EndTime: 2021-01-15 16:55:50
% DurationCPUTime: 0.53s
% Computational Cost: add. (628->113), mult. (1513->188), div. (0->0), fcn. (914->6), ass. (0->86)
t35 = sin(pkin(7)) * pkin(1) + qJ(3);
t44 = cos(pkin(8));
t46 = sin(qJ(4));
t47 = cos(qJ(4));
t42 = sin(pkin(8));
t24 = -cos(pkin(7)) * pkin(1) - t44 * pkin(3) - t42 * pkin(6) - pkin(2);
t88 = qJ(5) * t42;
t61 = -t24 + t88;
t99 = -t47 * t44 * t35 + t61 * t46;
t81 = qJD(4) * t47;
t65 = qJ(5) * t81;
t18 = t24 * qJD(1) + qJD(3);
t29 = t35 * qJD(1);
t21 = t42 * qJD(2) + t44 * t29;
t77 = qJD(1) * qJD(3);
t64 = t44 * t77;
t82 = qJD(4) * t46;
t73 = -t18 * t81 + t21 * t82 - t47 * t64;
t80 = qJD(5) * t46;
t86 = qJD(1) * t42;
t1 = (-t65 - t80) * t86 - t73;
t69 = t42 * t82;
t58 = qJD(1) * t69;
t28 = qJ(5) * t58;
t84 = qJD(3) * t46;
t67 = t44 * t84;
t92 = t42 * t47;
t50 = -qJD(5) * t92 - t67;
t53 = -t46 * t18 - t47 * t21;
t51 = t53 * qJD(4);
t2 = t50 * qJD(1) + t28 + t51;
t79 = t44 * qJD(1);
t34 = -qJD(4) + t79;
t62 = t47 * t18 - t46 * t21;
t66 = qJ(5) * t86;
t6 = -t47 * t66 + t62;
t3 = -t34 * pkin(4) + t6;
t7 = -t46 * t66 - t53;
t54 = t3 * t46 - t47 * t7;
t98 = -t54 * qJD(4) + t1 * t46 + t2 * t47;
t39 = t44 ^ 2;
t38 = t42 ^ 2;
t97 = 0.2e1 * t38;
t96 = t3 - t6;
t95 = t34 * t44;
t94 = t35 * t46;
t48 = qJD(1) ^ 2;
t93 = t38 * t48;
t83 = qJD(3) * t47;
t91 = t24 * t81 + t44 * t83;
t76 = qJD(1) * qJD(4);
t63 = t47 * t76;
t57 = t42 * t63;
t23 = pkin(4) * t57 + t42 * t77;
t90 = t38 + t39;
t40 = t46 ^ 2;
t41 = t47 ^ 2;
t89 = t40 - t41;
t87 = qJD(1) * t38;
t85 = qJD(1) * t46;
t37 = t44 * qJD(2);
t12 = qJD(5) - t37 + (pkin(4) * t85 + t29) * t42;
t78 = qJD(5) + t12;
t75 = t46 * t93;
t72 = t42 * t85;
t71 = t47 * t86;
t70 = t35 * t82;
t68 = t42 * t81;
t60 = qJD(1) * t90;
t59 = t34 * t69;
t55 = t3 * t47 + t46 * t7;
t20 = t42 * t29 - t37;
t52 = t20 * t42 + t21 * t44;
t49 = -t34 ^ 2 - t93;
t27 = t44 * t58;
t26 = (pkin(4) * t81 + qJD(3)) * t42;
t25 = t34 * t71;
t22 = (pkin(4) * t46 + t35) * t42;
t14 = (t34 - t79) * t68;
t13 = t27 - t59;
t11 = t49 * t47;
t10 = t49 * t46;
t8 = -t61 * t47 + (-pkin(4) - t94) * t44;
t5 = t99 * qJD(4) + t50;
t4 = -t42 * t80 + (-t44 * t94 - t47 * t88) * qJD(4) + t91;
t9 = [0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t60, (t35 * t60 + t52) * qJD(3), -0.2e1 * t38 * t46 * t63, t89 * t76 * t97, t27 + t59, (t34 + t79) * t68, 0, (t95 + (t97 + t39) * qJD(1)) * t84 + ((t18 * t44 + t24 * t34) * t46 + ((t87 + t95) * t35 + t52) * t47) * qJD(4), (-t44 * t70 + t91) * t34 - t73 * t44 - t20 * t69 + (-t70 + 0.2e1 * t83) * t87, -t2 * t44 - t5 * t34 + (t12 * t81 + t23 * t46 + (t22 * t81 + t26 * t46) * qJD(1)) * t42, t1 * t44 + t4 * t34 + (-t12 * t82 + t23 * t47 + (-t22 * t82 + t26 * t47) * qJD(1)) * t42, ((-t4 * t46 - t47 * t5 + (t46 * t8 + t47 * t99) * qJD(4)) * qJD(1) - t98) * t42, -t1 * t99 + t12 * t26 + t2 * t8 + t23 * t22 + t3 * t5 + t7 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, t14, t13, 0, -t23 * t44 + (-qJD(4) * t55 + t1 * t47 - t2 * t46) * t42; 0, 0, 0, 0, 0, -t90 * t48, -t52 * qJD(1), 0, 0, 0, 0, 0, t10, t11, t10, t11, 0, (-t12 * t42 + t44 * t54) * qJD(1) + t98; 0, 0, 0, 0, 0, 0, 0, t47 * t75, -t89 * t93, (-qJD(4) - t34) * t72, -t25 - t57, 0, t53 * t34 + t51 + (-t20 * t92 - t67) * qJD(1), t20 * t72 - t34 * t62 + t73, -t7 * t34 + t28 + (-qJD(4) * t18 - t64) * t46 + (-pkin(4) * t75 - qJD(4) * t21 - t78 * t86) * t47, -t41 * pkin(4) * t93 - t6 * t34 + (t46 * t78 + t65) * t86 + t73, (pkin(4) * qJD(4) - t96) * t72, t96 * t7 + (-t12 * t71 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 + t57, (-qJD(4) + t34) * t72, (-t40 - t41) * t93, t55 * t86 + t23;];
tauc_reg = t9;
