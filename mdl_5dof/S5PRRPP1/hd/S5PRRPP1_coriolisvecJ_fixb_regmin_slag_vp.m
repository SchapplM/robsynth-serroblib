% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:53
% EndTime: 2019-12-05 16:06:55
% DurationCPUTime: 0.41s
% Computational Cost: add. (743->126), mult. (1885->179), div. (0->0), fcn. (1195->4), ass. (0->77)
t82 = (qJD(2) * qJD(3));
t96 = -2 * t82;
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t62 = sin(qJ(3));
t63 = cos(qJ(3));
t46 = t60 * t63 + t61 * t62;
t42 = t46 * qJD(2);
t38 = t42 ^ 2;
t91 = t61 * t63;
t79 = qJD(2) * t91;
t85 = qJD(2) * t62;
t39 = t60 * t85 - t79;
t95 = -t39 ^ 2 - t38;
t87 = -qJ(4) - pkin(6);
t49 = t87 * t63;
t77 = t87 * t62;
t25 = -t60 * t49 - t61 * t77;
t74 = qJD(3) * t87;
t35 = t63 * qJD(4) + t62 * t74;
t81 = qJD(3) * qJD(1);
t24 = t35 * qJD(2) + t63 * t81;
t68 = -t62 * qJD(4) + t63 * t74;
t66 = t68 * qJD(2) - t62 * t81;
t3 = t60 * t24 - t61 * t66;
t94 = t3 * t25;
t45 = t60 * t62 - t91;
t93 = t3 * t45;
t36 = t62 * qJD(1) - qJD(2) * t49;
t92 = t60 * t36;
t28 = t61 * t36;
t65 = qJD(2) ^ 2;
t90 = t63 * t65;
t64 = qJD(3) ^ 2;
t89 = t64 * t62;
t88 = t64 * t63;
t4 = t61 * t24 + t60 * t66;
t34 = t63 * qJD(1) + qJD(2) * t77;
t31 = qJD(3) * pkin(3) + t34;
t12 = t60 * t31 + t28;
t86 = t62 ^ 2 - t63 ^ 2;
t84 = t62 * qJD(3);
t16 = t61 * t34 - t92;
t83 = qJD(5) - t16;
t80 = pkin(3) * t84;
t78 = -t63 * pkin(3) - pkin(2);
t76 = t62 * t82;
t75 = t63 * t82;
t73 = pkin(2) * t96;
t72 = t78 * qJD(2);
t11 = t61 * t31 - t92;
t41 = t46 * qJD(3);
t32 = qJD(2) * t41;
t50 = t60 * t76;
t33 = t61 * t75 - t50;
t53 = pkin(3) * t76;
t71 = t32 * pkin(4) - t33 * qJ(5) + t53;
t48 = qJD(4) + t72;
t10 = t39 * pkin(4) - t42 * qJ(5) + t48;
t70 = t10 * t42 + t3;
t44 = qJD(3) * t91 - t60 * t84;
t69 = -t46 * t32 + t45 * t33 - t44 * t39 + t41 * t42;
t15 = t60 * t35 - t61 * t68;
t17 = t61 * t35 + t60 * t68;
t26 = -t61 * t49 + t60 * t77;
t67 = t15 * t42 - t17 * t39 + t25 * t33 - t26 * t32 + t3 * t46;
t56 = -t61 * pkin(3) - pkin(4);
t54 = t60 * pkin(3) + qJ(5);
t18 = t45 * pkin(4) - t46 * qJ(5) + t78;
t14 = t60 * t34 + t28;
t13 = pkin(3) * t85 + t42 * pkin(4) + t39 * qJ(5);
t8 = qJD(3) * qJ(5) + t12;
t7 = -qJD(3) * pkin(4) + qJD(5) - t11;
t6 = t41 * pkin(4) - t44 * qJ(5) - t46 * qJD(5) + t80;
t2 = qJD(3) * qJD(5) + t4;
t1 = -t42 * qJD(5) + t71;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, t69, -t11 * t41 + t12 * t44 + t4 * t46 + t93, -t41 * qJD(3), t69, t44 * qJD(3), t2 * t46 + t7 * t41 + t8 * t44 + t93; 0, 0, 0, 0, 0.2e1 * t62 * t75, t86 * t96, t88, -t89, 0, -pkin(6) * t88 + t62 * t73, pkin(6) * t89 + t63 * t73, -t11 * t44 - t12 * t41 - t4 * t45 + t67, -t11 * t15 + t12 * t17 + t94 + t4 * t26 + (t48 + t72) * t80, -t15 * qJD(3) + t1 * t45 + t10 * t41 + t18 * t32 + t6 * t39, -t2 * t45 - t8 * t41 + t7 * t44 + t67, t17 * qJD(3) - t1 * t46 - t10 * t44 - t18 * t33 - t6 * t42, t1 * t18 + t10 * t6 + t7 * t15 + t8 * t17 + t2 * t26 + t94; 0, 0, 0, 0, -t62 * t90, t86 * t65, 0, 0, 0, t65 * pkin(2) * t62, pkin(2) * t90, (t12 - t14) * t42 + (-t11 + t16) * t39 + (-t32 * t60 - t33 * t61) * pkin(3), t11 * t14 - t12 * t16 + (-t3 * t61 + t4 * t60 - t48 * t85) * pkin(3), t14 * qJD(3) - t13 * t39 - t70, -t54 * t32 + t56 * t33 + (-t14 + t8) * t42 + (t7 - t83) * t39, -t10 * t39 + t13 * t42 + (0.2e1 * qJD(5) - t16) * qJD(3) + t4, -t10 * t13 - t7 * t14 + t2 * t54 + t3 * t56 + t83 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t11 * t42 + t12 * t39 + t53, 0.2e1 * t42 * qJD(3), t95, t50 + (t39 - t79) * qJD(3), t8 * t39 + (-qJD(5) - t7) * t42 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t39, -t50 + (t39 + t79) * qJD(3), -t38 - t64, -t8 * qJD(3) + t70;];
tauc_reg = t5;
