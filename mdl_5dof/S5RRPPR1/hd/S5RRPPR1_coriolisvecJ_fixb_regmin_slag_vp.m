% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:27
% EndTime: 2019-12-05 18:18:29
% DurationCPUTime: 0.40s
% Computational Cost: add. (570->94), mult. (1234->140), div. (0->0), fcn. (793->8), ass. (0->76)
t65 = sin(pkin(9));
t67 = cos(pkin(9));
t89 = t65 ^ 2 + t67 ^ 2;
t66 = sin(pkin(8));
t70 = sin(qJ(2));
t88 = pkin(1) * qJD(1);
t84 = t70 * t88;
t51 = t66 * t84;
t68 = cos(pkin(8));
t72 = cos(qJ(2));
t83 = t72 * t88;
t37 = t68 * t83 - t51;
t27 = t37 * qJD(2);
t64 = qJD(1) + qJD(2);
t19 = t64 * qJD(4) + t27;
t102 = t19 * t89;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t44 = t71 * t65 + t69 * t67;
t30 = t44 * t64;
t101 = t64 * t89;
t100 = t37 - qJD(4);
t99 = t67 * pkin(4);
t98 = t68 * pkin(2);
t87 = pkin(1) * qJD(2);
t94 = t68 * t70;
t36 = (t66 * t72 + t94) * t87;
t26 = qJD(1) * t36;
t39 = t44 * qJD(5);
t92 = t71 * t67;
t93 = t69 * t65;
t43 = -t92 + t93;
t45 = t64 * pkin(2) + t83;
t20 = t68 * t45 - t51;
t76 = qJD(4) - t20;
t82 = -pkin(3) - t99;
t9 = t82 * t64 + t76;
t97 = t26 * t43 + t9 * t39;
t38 = t43 * qJD(5);
t96 = t26 * t44 - t9 * t38;
t95 = t64 * t65;
t52 = t68 * t84;
t21 = t66 * t45 + t52;
t57 = t72 * pkin(1) + pkin(2);
t90 = pkin(1) * t94 + t66 * t57;
t54 = t66 * t70 * pkin(1);
t86 = t64 * t93;
t85 = t64 * t92;
t80 = t68 * t57 - t54;
t78 = -pkin(3) - t80;
t77 = t89 * (t64 * qJ(4) + t21);
t75 = (-qJD(2) + t64) * t88;
t74 = (-qJD(1) - t64) * t87;
t73 = t68 * t72 * t87 - qJD(2) * t54;
t60 = t67 * pkin(7);
t56 = t66 * pkin(2) + qJ(4);
t47 = t82 - t98;
t46 = qJD(5) * t85;
t41 = t67 * t56 + t60;
t40 = (-pkin(7) - t56) * t65;
t35 = t66 * t83 + t52;
t34 = t39 * qJD(5);
t33 = t38 * qJD(5);
t32 = qJ(4) + t90;
t31 = qJD(4) + t73;
t28 = -t85 + t86;
t25 = t78 - t99;
t24 = t64 * t39;
t23 = -qJD(5) * t86 + t46;
t22 = t26 * t65;
t18 = t67 * t32 + t60;
t17 = (-pkin(7) - t32) * t65;
t13 = -t64 * pkin(3) + t76;
t2 = t23 * t44 - t30 * t38;
t1 = -t23 * t43 - t44 * t24 + t38 * t28 - t30 * t39;
t3 = [0, 0, 0, 0, t70 * t74, t72 * t74, -t20 * t36 + t21 * t73 - t26 * t80 + t27 * t90, (-t36 * t64 - t26) * t67, t36 * t95 + t22, t31 * t101 + t102, t102 * t32 + t13 * t36 + t26 * t78 + t77 * t31, t2, t1, -t33, -t34, 0, t25 * t24 + t36 * t28 + ((-t17 * t69 - t18 * t71) * qJD(5) - t44 * t31) * qJD(5) + t97, t25 * t23 + t36 * t30 + ((-t17 * t71 + t18 * t69) * qJD(5) + t43 * t31) * qJD(5) + t96; 0, 0, 0, 0, t70 * t75, t72 * t75, t20 * t35 - t21 * t37 + (-t26 * t68 + t27 * t66) * pkin(2), (t35 * t64 - t26) * t67, -t35 * t95 + t22, -t100 * t101 + t102, t26 * (-pkin(3) - t98) - t13 * t35 + t56 * t102 - t100 * t77, t2, t1, -t33, -t34, 0, t47 * t24 - t35 * t28 + ((-t40 * t69 - t41 * t71) * qJD(5) + t100 * t44) * qJD(5) + t97, t47 * t23 - t35 * t30 + ((-t40 * t71 + t41 * t69) * qJD(5) - t100 * t43) * qJD(5) + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * t64 ^ 2, -t77 * t64 + t26, 0, 0, 0, 0, 0, 0.2e1 * t30 * qJD(5), t46 + (-t28 - t86) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t28, -t28 ^ 2 + t30 ^ 2, t46 + (t28 - t86) * qJD(5), 0, 0, -t44 * t19 - t9 * t30, t43 * t19 + t9 * t28;];
tauc_reg = t3;
