% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:31
% DurationCPUTime: 0.75s
% Computational Cost: add. (542->160), mult. (1332->253), div. (0->0), fcn. (905->8), ass. (0->103)
t37 = sin(qJ(4));
t78 = t37 * qJD(2);
t29 = qJD(5) + t78;
t108 = -qJD(5) + t29;
t39 = cos(qJ(5));
t40 = cos(qJ(4));
t85 = qJD(2) * t40;
t66 = t39 * t85;
t36 = sin(qJ(5));
t79 = t36 * qJD(4);
t24 = t66 + t79;
t75 = qJD(2) * qJD(4);
t97 = t36 * t37;
t14 = t24 * qJD(5) - t75 * t97;
t42 = -pkin(2) - pkin(7);
t41 = cos(qJ(2));
t34 = sin(pkin(5));
t89 = qJD(1) * t34;
t67 = t41 * t89;
t56 = qJD(3) - t67;
t18 = t42 * qJD(2) + t56;
t35 = cos(pkin(5));
t38 = sin(qJ(2));
t86 = qJD(2) * t34;
t72 = t38 * t86;
t49 = -qJD(4) * t35 + t72;
t84 = qJD(4) * t37;
t87 = qJD(1) * t40;
t4 = t18 * t84 - t49 * t87;
t107 = t4 * t36;
t106 = t4 * t39;
t81 = qJD(5) * t36;
t70 = t40 * t81;
t77 = t39 * qJD(4);
t48 = -t37 * t77 - t70;
t13 = t48 * qJD(2) + qJD(5) * t77;
t105 = t13 * t36;
t22 = t36 * t85 - t77;
t104 = t22 * t29;
t103 = t24 * t29;
t68 = t38 * t89;
t76 = qJD(2) * qJ(3);
t26 = t68 + t76;
t102 = t26 * t41;
t101 = t29 * t36;
t100 = t34 * t38;
t99 = t34 * t41;
t44 = qJD(2) ^ 2;
t98 = t34 * t44;
t96 = t36 * t42;
t95 = t39 * t29;
t94 = t39 * t42;
t93 = t40 * t13;
t33 = t40 ^ 2;
t92 = t37 ^ 2 - t33;
t43 = qJD(4) ^ 2;
t91 = -t43 - t44;
t90 = qJD(2) * pkin(2);
t88 = qJD(1) * t37;
t83 = qJD(4) * t40;
t82 = qJD(4) * t42;
t80 = qJD(5) * t39;
t74 = t38 * t98;
t73 = t41 * t98;
t71 = t41 * t86;
t69 = t29 * t80;
t11 = t37 * t18 + t35 * t87;
t7 = qJD(4) * pkin(8) + t11;
t65 = t29 * t42 + t7;
t64 = t40 * t75;
t63 = t29 + t78;
t62 = t40 * t72;
t61 = t37 * t72;
t60 = qJD(5) * t37 + qJD(2);
t59 = pkin(4) * t40 + pkin(8) * t37;
t19 = t59 * qJD(4) + qJD(3);
t58 = -t19 + t67;
t57 = -t26 + t68;
t27 = t37 * pkin(4) - t40 * pkin(8) + qJ(3);
t15 = t27 * qJD(2) + t68;
t1 = t39 * t15 - t36 * t7;
t2 = t36 * t15 + t39 * t7;
t55 = qJD(2) * t33 - t29 * t37;
t51 = -t40 * t18 + t35 * t88;
t6 = -qJD(4) * pkin(4) + t51;
t54 = -pkin(8) * t83 + t37 * t6;
t17 = t35 * t40 - t37 * t99;
t53 = t39 * t100 - t17 * t36;
t52 = t36 * t100 + t17 * t39;
t16 = t35 * t37 + t40 * t99;
t50 = t57 * qJD(2);
t47 = t57 - t76;
t20 = (qJD(3) + t67) * qJD(2);
t46 = t56 * qJD(2) - t42 * t43 + t20;
t3 = t18 * t83 + t49 * t88;
t45 = -qJD(4) * t6 - qJD(5) * t15 + t29 * t68 - t3;
t25 = t59 * qJD(2);
t21 = t56 - t90;
t12 = (t19 + t67) * qJD(2);
t9 = t17 * qJD(4) - t62;
t8 = -t16 * qJD(4) + t61;
t5 = t39 * t12;
t10 = [0, 0, -t74, -t73, t74, t73, (t20 * t38 + (t102 + (t21 - t67) * t38) * qJD(2)) * t34, 0, 0, 0, 0, 0, t37 * t73 + (-t9 + t62) * qJD(4), t40 * t73 + (-t8 - t61) * qJD(4), 0, 0, 0, 0, 0, (-t52 * qJD(5) - t8 * t36 + t39 * t71) * t29 + t53 * t64 + t9 * t22 + t16 * t14, -(t53 * qJD(5) + t36 * t71 + t8 * t39) * t29 - t52 * t64 + t9 * t24 + t16 * t13; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t20 * qJ(3) + t26 * qJD(3) + (-t102 + (-t21 - t90) * t38) * t89, -0.2e1 * t37 * t64, 0.2e1 * t92 * t75, -t43 * t37, -t43 * t40, 0, t46 * t37 - t47 * t83, t46 * t40 + t47 * t84, t48 * t24 + t39 * t93, (t22 * t39 + t24 * t36) * t84 + (-t105 - t14 * t39 + (t22 * t36 - t24 * t39) * qJD(5)) * t40, -t29 * t70 + t13 * t37 + (t24 * t40 + t55 * t39) * qJD(4), -t40 * t69 - t14 * t37 + (-t22 * t40 - t55 * t36) * qJD(4), t63 * t83, (-t27 * t81 - t58 * t39) * t29 + (t22 * t82 + t45 * t36 - t65 * t80 + t5) * t37 + (t22 * t68 + t6 * t80 - t42 * t14 + t107 + (-t29 * t96 + (t39 * t27 - t37 * t96) * qJD(2) + t1) * qJD(4)) * t40, (-t27 * t80 + t58 * t36) * t29 + (t24 * t82 + (t65 * qJD(5) - t12) * t36 + t45 * t39) * t37 + (t24 * t68 - t6 * t81 - t42 * t13 + t106 + (-t29 * t94 - (t36 * t27 + t37 * t94) * qJD(2) - t2) * qJD(4)) * t40; 0, 0, 0, 0, 0, -t44, t50, 0, 0, 0, 0, 0, t91 * t37, t91 * t40, 0, 0, 0, 0, 0, -t40 * t14 - t60 * t95 + (-t63 * t40 * t36 + t22 * t37) * qJD(4), -t93 + t60 * t101 + (-t40 * t95 + (t24 - t66) * t37) * qJD(4); 0, 0, 0, 0, 0, 0, 0, t40 * t44 * t37, -t92 * t44, 0, 0, 0, t40 * t50, -t57 * t78, t24 * t95 + t105, (t13 - t104) * t39 + (-t14 - t103) * t36, t69 + (t37 * t95 + (-t24 + t79) * t40) * qJD(2), -t29 * t81 + (-t29 * t97 + (t22 + t77) * t40) * qJD(2), -t29 * t85, -pkin(4) * t14 - t106 - (t39 * t25 + t36 * t51) * t29 - t11 * t22 + (-pkin(8) * t95 + t6 * t36) * qJD(5) + (-t1 * t40 + t54 * t36) * qJD(2), -pkin(4) * t13 + t107 + (t36 * t25 - t39 * t51) * t29 - t11 * t24 + (pkin(8) * t101 + t6 * t39) * qJD(5) + (t2 * t40 + t54 * t39) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t22, -t22 ^ 2 + t24 ^ 2, t13 + t104, t103 - t14, t64, t108 * t2 - t6 * t24 - t36 * t3 + t5, t108 * t1 - t36 * t12 + t6 * t22 - t39 * t3;];
tauc_reg = t10;
