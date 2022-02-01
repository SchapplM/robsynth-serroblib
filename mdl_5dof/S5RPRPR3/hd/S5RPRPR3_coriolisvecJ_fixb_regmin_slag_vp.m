% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:56
% DurationCPUTime: 0.40s
% Computational Cost: add. (588->87), mult. (1258->146), div. (0->0), fcn. (699->8), ass. (0->77)
t41 = qJD(1) + qJD(3);
t36 = cos(pkin(8)) * pkin(1) + pkin(2);
t32 = t36 * qJD(1);
t49 = sin(qJ(3));
t95 = pkin(1) * sin(pkin(8));
t69 = qJD(3) * t95;
t62 = qJD(1) * t69;
t51 = cos(qJ(3));
t76 = qJD(3) * t51;
t60 = t32 * t76 - t49 * t62;
t13 = t41 * qJD(4) + t60;
t44 = sin(pkin(9));
t39 = t44 ^ 2;
t46 = cos(pkin(9));
t78 = t46 ^ 2 + t39;
t100 = t78 * t13;
t81 = t49 * t32;
t19 = qJD(3) * t81 + t51 * t62;
t54 = t49 * t36 + t51 * t95;
t23 = t54 * qJD(3);
t99 = -t23 * t41 - t19;
t70 = qJD(1) * t95;
t21 = t51 * t70 + t81;
t98 = t21 * t41 - t19;
t97 = -t51 * t36 + t49 * t95;
t20 = t51 * t32 - t49 * t70;
t57 = qJD(4) - t20;
t29 = -t46 * pkin(4) - t44 * pkin(7) - pkin(3);
t3 = t29 * t41 + t57;
t48 = sin(qJ(5));
t50 = cos(qJ(5));
t16 = t41 * qJ(4) + t21;
t6 = t44 * qJD(2) + t46 * t16;
t82 = t46 * t50;
t86 = t39 * t50;
t96 = (t13 * t82 + t48 * t19 + (t50 * t3 - t48 * t6) * qJD(5)) * t46 + t13 * t86;
t5 = -t46 * qJD(2) + t44 * t16;
t94 = t5 * t44;
t55 = t36 * t76 - t49 * t69;
t22 = qJD(4) + t55;
t91 = t22 * t41;
t38 = t41 ^ 2;
t89 = t39 * t38;
t88 = t39 * t41;
t87 = t39 * t48;
t85 = t41 * t44;
t84 = t46 * t41;
t83 = t46 * t48;
t33 = -qJD(5) + t84;
t80 = t50 * t33;
t77 = t48 ^ 2 - t50 ^ 2;
t75 = qJD(5) * t48;
t74 = qJD(5) * t50;
t73 = qJD(5) + t33;
t72 = t5 * t85;
t71 = t41 * t86;
t68 = qJD(5) * t88;
t67 = t44 * t75;
t66 = t44 * t74;
t59 = -t48 * t3 - t50 * t6;
t2 = t59 * qJD(5) - t13 * t83 + t50 * t19;
t64 = t13 * t87 - t2 * t46 + t5 * t66;
t63 = t33 * t67;
t61 = t73 * t85;
t58 = t6 * t46 + t94;
t56 = t33 * t46 + t88;
t53 = qJ(4) * t74 + t57 * t48;
t52 = -t33 ^ 2 - t89;
t27 = t67 * t84;
t25 = -0.2e1 * t50 * t48 * t68;
t24 = qJ(4) + t54;
t18 = 0.2e1 * t77 * t68;
t17 = t29 + t97;
t14 = -t41 * pkin(3) + t57;
t12 = (t33 + t84) * t66;
t11 = t27 + t63;
t1 = [0, 0, 0, 0, 0, t99, -t55 * t41 - t60, t99 * t46, t78 * t91 + t100, t19 * (-pkin(3) + t97) + t14 * t23 + t58 * t22 + t24 * t100, t25, t18, t11, t12, 0, -(-t22 * t83 + t50 * t23) * t33 + t87 * t91 + (-(-t17 * t48 - t24 * t82) * t33 + t24 * t71) * qJD(5) + t64, (t22 * t82 + t48 * t23) * t33 + t22 * t71 + (t17 * t80 + (-t56 * t24 - t94) * t48) * qJD(5) + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t33 - t84) * t66, t27 - t63; 0, 0, 0, 0, 0, t98, t20 * t41 - t60, t98 * t46, t57 * t41 * t78 + t100, -t19 * pkin(3) + qJ(4) * t100 - t14 * t21 + t57 * t58, t25, t18, t11, t12, 0, (t50 * t21 + t29 * t75 + t53 * t46) * t33 + t53 * t88 + t64, -t48 * t21 * t33 + t57 * t56 * t50 + (t29 * t80 + (-t56 * qJ(4) - t94) * t48) * qJD(5) + t96; 0, 0, 0, 0, 0, 0, 0, 0, -t78 * t38, -t58 * t41 + t19, 0, 0, 0, 0, 0, t52 * t48, t52 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t48 * t86, -t77 * t89, -t48 * t61, -t50 * t61, 0, t59 * t33 - t50 * t72 + t2, (-t46 * t13 - t73 * t3) * t50 + (t73 * t6 - t19 + t72) * t48;];
tauc_reg = t1;
