% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:50
% DurationCPUTime: 0.54s
% Computational Cost: add. (1069->113), mult. (1617->160), div. (0->0), fcn. (806->6), ass. (0->88)
t54 = sin(pkin(8));
t55 = cos(pkin(8));
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t90 = pkin(1) * qJD(1);
t23 = (t54 * t57 + t55 * t59) * t90;
t104 = -t55 * qJD(3) + t23;
t50 = qJD(1) + qJD(2);
t107 = t104 * t50;
t60 = -pkin(2) - pkin(3);
t93 = t55 * qJ(3) + t54 * t60;
t31 = -pkin(7) + t93;
t61 = qJD(5) ^ 2;
t81 = t59 * t90;
t82 = t57 * t90;
t22 = t54 * t81 - t55 * t82;
t66 = (t54 * qJD(3) - t22) * t50;
t106 = -t31 * t61 + t66;
t89 = pkin(1) * qJD(2);
t79 = qJD(1) * t89;
t42 = t59 * t79;
t46 = t50 * qJD(3);
t28 = t42 + t46;
t74 = t57 * t79;
t14 = t55 * t28 + t54 * t74;
t58 = cos(qJ(5));
t56 = sin(qJ(5));
t71 = qJD(3) - t81;
t21 = t50 * t60 + t71;
t32 = t50 * qJ(3) + t82;
t11 = t54 * t21 + t55 * t32;
t9 = -t50 * pkin(7) + t11;
t6 = t58 * qJD(4) - t56 * t9;
t2 = t6 * qJD(5) + t58 * t14;
t7 = t56 * qJD(4) + t58 * t9;
t3 = -t7 * qJD(5) - t56 * t14;
t105 = t2 * t58 - t3 * t56;
t49 = t50 ^ 2;
t10 = t55 * t21 - t54 * t32;
t8 = t50 * pkin(4) - t10;
t102 = t50 * t8;
t101 = t56 * t7;
t83 = t59 * t89;
t41 = qJD(3) + t83;
t84 = t57 * t89;
t18 = t54 * t41 - t55 * t84;
t100 = t18 * t50;
t19 = t55 * t41 + t54 * t84;
t99 = t19 * t50;
t97 = t50 * t56;
t96 = t55 * t49;
t48 = t61 * t56;
t95 = t61 * t58;
t80 = -t59 * pkin(1) - pkin(2);
t43 = -pkin(3) + t80;
t44 = t57 * pkin(1) + qJ(3);
t94 = t54 * t43 + t55 * t44;
t52 = t56 ^ 2;
t53 = t58 ^ 2;
t92 = t52 - t53;
t91 = t52 + t53;
t88 = qJD(5) * t58;
t86 = t56 * t49 * t58;
t85 = 0.2e1 * qJD(5) * t55;
t78 = -t14 + t102;
t13 = t54 * t28 - t55 * t74;
t75 = t88 * t97;
t73 = t56 * t6 - t58 * t7;
t72 = qJD(5) * t101 + t6 * t88 - t105;
t70 = t10 * t54 - t11 * t55;
t16 = -pkin(7) + t94;
t69 = -t16 * t61 + t100;
t68 = t55 * t43 - t54 * t44;
t67 = -t54 * qJ(3) + t55 * t60;
t65 = t50 * t81 - t42;
t15 = pkin(4) - t68;
t64 = qJD(5) * (-t15 * t50 - t19 - t8);
t30 = pkin(4) - t67;
t63 = qJD(5) * (-t30 * t50 + t104 - t8);
t62 = (-t58 * t6 - t101) * qJD(5) + t105;
t34 = 0.2e1 * t75;
t33 = -0.2e1 * t75;
t29 = -t50 * pkin(2) + t71;
t26 = (-qJD(1) - t50) * t84;
t25 = (-qJD(2) + t50) * t82;
t20 = -0.2e1 * t92 * t50 * qJD(5);
t12 = t13 * t58;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t50 * t83 - t42, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t41 * t50 + t28, t28 * t44 + t32 * t41 + (qJD(1) * t80 + t29) * t84, 0, 0, 0, 0, 0, 0, t13 + t100, t14 + t99, 0, -t10 * t18 + t11 * t19 - t13 * t68 + t14 * t94, t34, t20, -t95, t33, t48, 0, t56 * t64 + t58 * t69 + t12, (-t13 - t69) * t56 + t58 * t64, -t91 * t99 + t72, t13 * t15 + t16 * t62 + t8 * t18 - t19 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t65, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0.2e1 * t46 - t65, t28 * qJ(3) + t32 * qJD(3) + (-t32 * t59 + (-pkin(2) * qJD(2) - t29) * t57) * t90, 0, 0, 0, 0, 0, 0, t13 + t66, t14 - t107, 0, -qJD(3) * t70 + t10 * t22 - t11 * t23 - t13 * t67 + t14 * t93, t34, t20, -t95, t33, t48, 0, t106 * t58 + t56 * t63 + t12, (-t13 - t106) * t56 + t58 * t63, t91 * t107 + t72, t13 * t30 - t8 * t22 + t73 * t23 + (t54 * t8 - t55 * t73) * qJD(3) + t62 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t32 * t50 + t74, 0, 0, 0, 0, 0, 0, -t54 * t49, -t96, 0, -t13 * t55 + t14 * t54 + t50 * t70, 0, 0, 0, 0, 0, 0, -t54 * t95 + (-t54 * t58 * t50 + t56 * t85) * t50, t54 * t48 + (t54 * t97 + t58 * t85) * t50, t91 * t96, (t50 * t73 - t13) * t55 + (t62 - t102) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t95, 0, -qJD(5) * t73 + t2 * t56 + t3 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t92 * t49, 0, t86, 0, 0, t78 * t56, t78 * t58, 0, 0;];
tauc_reg = t1;
