% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP2
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
% tauc_reg [4x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:29
% EndTime: 2021-01-15 11:04:30
% DurationCPUTime: 0.39s
% Computational Cost: add. (469->106), mult. (860->158), div. (0->0), fcn. (399->4), ass. (0->83)
t82 = -qJ(4) - pkin(6);
t44 = sin(qJ(3));
t41 = qJD(1) + qJD(2);
t45 = sin(qJ(2));
t77 = pkin(1) * qJD(1);
t67 = t45 * t77;
t61 = -t82 * t41 + t67;
t9 = t61 * t44;
t46 = cos(qJ(3));
t10 = t61 * t46;
t75 = qJD(3) * pkin(3);
t6 = -t9 + t75;
t91 = t6 + t9;
t42 = t44 ^ 2;
t90 = pkin(3) * t42;
t89 = t46 * pkin(3);
t47 = cos(qJ(2));
t88 = t47 * pkin(1);
t36 = -pkin(2) - t89;
t66 = t47 * t77;
t11 = t36 * t41 + qJD(4) - t66;
t76 = pkin(1) * qJD(2);
t62 = qJD(1) * t76;
t33 = t45 * t62;
t73 = qJD(3) * t44;
t65 = t41 * t73;
t16 = pkin(3) * t65 + t33;
t72 = qJD(3) * t46;
t87 = t11 * t72 + t16 * t44;
t40 = t41 ^ 2;
t86 = t40 * t46;
t85 = t41 * t44;
t84 = t41 * t46;
t48 = qJD(3) ^ 2;
t83 = t48 * t44;
t38 = t48 * t46;
t25 = -t41 * pkin(2) - t66;
t81 = t25 * t72 + t44 * t33;
t56 = qJD(3) * t66;
t57 = t41 * t67;
t80 = t44 * t56 + t46 * t57;
t43 = t46 ^ 2;
t79 = -t42 - t43;
t78 = t42 - t43;
t34 = t45 * pkin(1) + pkin(6);
t74 = -qJ(4) - t34;
t71 = qJD(3) * t47;
t70 = -qJD(1) - t41;
t69 = t47 * t76;
t68 = t45 * t76;
t7 = t11 * t73;
t64 = t41 * t72;
t55 = t47 * t62;
t50 = qJD(4) * t41 + t55;
t52 = qJD(3) * t61;
t2 = -t44 * t52 + t50 * t46;
t3 = -t50 * t44 - t46 * t52;
t63 = t2 * t46 - t3 * t44;
t60 = qJD(3) * t82;
t59 = 0.2e1 * t64;
t58 = qJD(3) * t74;
t54 = -t10 * t46 + t44 * t6;
t53 = -t10 * t44 - t46 * t6;
t51 = -t25 * t41 - t55;
t49 = -t55 + (-qJD(4) - t11) * t41;
t39 = t46 * qJ(4);
t37 = t46 * qJD(4);
t35 = -pkin(2) - t88;
t31 = t46 * pkin(6) + t39;
t30 = t82 * t44;
t28 = t46 * t56;
t26 = t36 - t88;
t23 = t44 * t59;
t21 = pkin(3) * t73 + t68;
t20 = t46 * t34 + t39;
t19 = t74 * t44;
t17 = t25 * t73;
t15 = -t44 * qJD(4) + t46 * t60;
t14 = t44 * t60 + t37;
t13 = -0.2e1 * t78 * t41 * qJD(3);
t5 = (-qJD(4) - t69) * t44 + t46 * t58;
t4 = t44 * t58 + t46 * t69 + t37;
t1 = [0, 0, 0, 0, -t41 * t68 - t33, t70 * t69, t23, t13, t38, -t83, 0, t35 * t65 - t34 * t38 + t17 + (t70 * t46 * t45 - t44 * t71) * t76, t35 * t64 + t34 * t83 + (t45 * t85 - t46 * t71) * t76 + t81, t7 + (-t21 * t41 - t16) * t46 + (t26 * t85 + t5) * qJD(3), t21 * t85 + (t26 * t84 - t4) * qJD(3) + t87, (t4 * t46 - t44 * t5) * t41 + ((-t19 * t46 - t20 * t44) * t41 + t53) * qJD(3) + t63, t10 * t4 + t11 * t21 + t16 * t26 + t3 * t19 + t2 * t20 + t6 * t5; 0, 0, 0, 0, -t33 + t57, (-qJD(2) + t41) * t66, t23, t13, t38, -t83, 0, -pkin(2) * t65 + t17 + (-pkin(6) * t48 - t33) * t46 + t80, pkin(6) * t83 + t28 + (-pkin(2) * t72 - t44 * t67) * t41 + t81, -t16 * t46 + t7 + (t15 + (t36 - t89) * t85) * qJD(3) + t80, -t44 * t57 + t28 + (-t14 + (t36 * t46 + t90) * t41) * qJD(3) + t87, t53 * qJD(3) + (t14 * t46 - t15 * t44 + (-t30 * t46 - t31 * t44) * qJD(3) + t79 * t66) * t41 + t63, pkin(3) * t7 + t10 * t14 + t6 * t15 + t16 * t36 + t2 * t31 + t3 * t30 + (-t11 * t45 + t54 * t47) * t77; 0, 0, 0, 0, 0, 0, -t44 * t86, t78 * t40, 0, 0, 0, t51 * t44, t51 * t46, (pkin(3) * t86 + t49) * t44, -t40 * t90 + t49 * t46, (-t75 + t91) * t84, t91 * t10 + (-t11 * t85 + t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65, t59, t79 * t40, t54 * t41 + t16;];
tauc_reg = t1;
