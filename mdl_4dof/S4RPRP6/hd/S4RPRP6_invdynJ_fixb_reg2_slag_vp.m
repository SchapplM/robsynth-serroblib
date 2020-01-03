% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:16
% EndTime: 2019-12-31 16:46:17
% DurationCPUTime: 0.54s
% Computational Cost: add. (397->138), mult. (744->154), div. (0->0), fcn. (345->4), ass. (0->92)
t54 = -pkin(1) - pkin(5);
t24 = t54 * qJD(1) + qJD(2);
t77 = qJ(4) * qJD(1);
t106 = t24 - t77;
t52 = cos(qJ(3));
t9 = t106 * t52;
t50 = sin(qJ(3));
t99 = t50 * pkin(3);
t31 = qJ(2) + t99;
t84 = qJD(1) * t31;
t16 = qJD(4) + t84;
t53 = cos(qJ(1));
t42 = g(2) * t53;
t51 = sin(qJ(1));
t43 = g(1) * t51;
t92 = t43 - t42;
t105 = -t16 * qJD(1) - t92;
t45 = qJDD(1) * qJ(2);
t46 = qJD(1) * qJD(2);
t73 = qJD(1) * qJD(3);
t67 = t52 * t73;
t75 = t50 * qJDD(1);
t4 = qJDD(4) + t45 + t46 + (t67 + t75) * pkin(3);
t61 = g(1) * t53 + g(2) * t51;
t104 = t4 - t61;
t23 = qJDD(1) * t54 + qJDD(2);
t47 = t50 ^ 2;
t48 = t52 ^ 2;
t90 = t47 + t48;
t69 = t90 * t23;
t71 = 0.2e1 * t46;
t103 = 0.2e1 * t45 + t71 - t61;
t74 = qJ(4) * qJDD(1);
t102 = (qJD(4) + t16) * qJD(1) + t74;
t88 = qJD(3) * pkin(3);
t5 = t9 + t88;
t100 = t5 - t9;
t40 = g(3) * t50;
t98 = t50 * t24;
t97 = t50 * t51;
t56 = qJD(1) ^ 2;
t96 = t56 * t50;
t95 = g(1) * t97 + g(3) * t52;
t94 = (t71 + t45) * qJ(2);
t93 = t53 * pkin(1) + t51 * qJ(2);
t91 = t47 - t48;
t55 = qJD(3) ^ 2;
t89 = -t55 - t56;
t87 = t56 * qJ(2);
t86 = qJ(4) - t54;
t85 = pkin(1) * qJDD(1);
t18 = t86 * t52;
t83 = qJD(3) * t18;
t82 = qJD(3) * t50;
t81 = qJD(3) * t52;
t80 = qJDD(3) * pkin(3);
t76 = qJDD(3) * t50;
t35 = t52 * qJDD(1);
t14 = t52 * t23;
t72 = t52 * t42 + t14 + t40;
t70 = -t23 - t42;
t68 = t50 * t73;
t66 = t16 + t84;
t19 = t90 * qJDD(1);
t64 = qJDD(2) - t85;
t63 = t50 * t67;
t62 = -t87 - t92;
t60 = -qJD(1) * qJD(4) - t74;
t59 = 0.2e1 * qJ(2) * t73 + qJDD(3) * t54;
t25 = pkin(3) * t81 + qJD(2);
t58 = qJD(1) * t25 + qJDD(1) * t31 + t104;
t57 = -t54 * t55 + t103;
t49 = -qJ(4) - pkin(5);
t38 = t53 * qJ(2);
t36 = qJDD(3) * t52;
t27 = t52 * t96;
t26 = qJ(4) * t68;
t22 = t91 * t56;
t21 = -t55 * t50 + t36;
t20 = -t55 * t52 - t76;
t17 = t86 * t50;
t13 = t48 * qJDD(1) - 0.2e1 * t63;
t12 = t47 * qJDD(1) + 0.2e1 * t63;
t11 = t89 * t50 + t36;
t10 = t89 * t52 - t76;
t8 = -t50 * t77 + t98;
t7 = -t50 * qJD(4) - t83;
t6 = -t52 * qJD(4) + t86 * t82;
t3 = -0.2e1 * t50 * t35 + 0.2e1 * t91 * t73;
t2 = t106 * t81 + (t23 + t60) * t50;
t1 = -t24 * t82 + t60 * t52 + t14 + t26 + t80;
t15 = [0, 0, 0, 0, 0, qJDD(1), t92, t61, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t85 - t92, t103, -t64 * pkin(1) - g(1) * (-t51 * pkin(1) + t38) - g(2) * t93 + t94, t13, t3, t21, t12, t20, 0, t57 * t50 + t59 * t52, -t59 * t50 + t57 * t52, -t54 * t19 - t69 + t92, -g(1) * (t54 * t51 + t38) - g(2) * (t53 * pkin(5) + t93) + t54 * t69 + t94, t13, t3, t21, t12, t20, 0, -t18 * qJDD(3) + (t66 * t52 + t6) * qJD(3) + t58 * t50, t17 * qJDD(3) + (-t66 * t50 - t7) * qJD(3) + t58 * t52, (-t8 * qJD(3) + qJDD(1) * t18 - t1 + (qJD(3) * t17 - t6) * qJD(1)) * t52 + (qJD(3) * t5 + qJDD(1) * t17 - t2 + (-t7 - t83) * qJD(1)) * t50 + t92, -t2 * t17 + t8 * t7 - t1 * t18 + t5 * t6 + t4 * t31 + t16 * t25 - g(1) * (t53 * t99 + t38 + (-pkin(1) + t49) * t51) - g(2) * (pkin(3) * t97 - t53 * t49 + t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t56, t62 + t64, 0, 0, 0, 0, 0, 0, t11, t10, -t19, t69 + t62, 0, 0, 0, 0, 0, 0, t11, t10, -t19, t1 * t52 + t2 * t50 + (-t5 * t50 + t52 * t8) * qJD(3) + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t22, t35, -t27, -t75, qJDD(3), (-t87 - t43) * t52 + t72, (t70 + t87) * t50 + t95, 0, 0, t27, -t22, t35, -t27, -t75, qJDD(3), 0.2e1 * t80 + t26 + (t8 - t98) * qJD(3) + (-pkin(3) * t96 - t102 - t43) * t52 + t72, -t48 * t56 * pkin(3) + (t102 + t70) * t50 + t95, -pkin(3) * t35 + (t88 - t100) * t50 * qJD(1), t100 * t8 + (t105 * t52 + t1 + t40) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67 + t75, t35 - 0.2e1 * t68, -t90 * t56, (t5 * t52 + t50 * t8) * qJD(1) + t104;];
tau_reg = t15;
