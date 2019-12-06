% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:45
% EndTime: 2019-12-05 16:00:48
% DurationCPUTime: 0.53s
% Computational Cost: add. (391->112), mult. (912->165), div. (0->0), fcn. (595->6), ass. (0->83)
t44 = sin(qJ(5));
t45 = sin(qJ(4));
t47 = cos(qJ(5));
t48 = cos(qJ(4));
t25 = t44 * t48 + t47 * t45;
t41 = qJD(4) + qJD(5);
t58 = t25 * t41;
t5 = qJD(2) * t58;
t96 = qJD(5) - t41;
t78 = qJD(4) * t48;
t98 = qJD(5) * t48 + t78;
t79 = qJD(4) * t45;
t97 = qJD(5) * t45 + t79;
t50 = -pkin(2) - pkin(6);
t95 = pkin(7) - t50;
t55 = t97 * t44;
t90 = t47 * t48;
t61 = t41 * t90;
t11 = t61 - t55;
t94 = t11 * t41;
t80 = qJD(2) * t48;
t69 = t47 * t80;
t81 = qJD(2) * t45;
t70 = t44 * t81;
t18 = -t69 + t70;
t19 = t25 * qJD(2);
t93 = t18 * t19;
t73 = qJD(2) * qJ(3);
t46 = sin(qJ(2));
t75 = t46 * qJD(1);
t31 = t73 + t75;
t49 = cos(qJ(2));
t92 = t31 * t49;
t74 = t49 * qJD(1);
t62 = qJD(3) - t74;
t24 = t50 * qJD(2) + t62;
t14 = -pkin(7) * t81 + t45 * t24;
t91 = t47 * t14;
t52 = qJD(2) ^ 2;
t89 = t52 * t46;
t88 = t52 * t49;
t72 = qJD(2) * qJD(4);
t65 = t45 * t72;
t87 = -qJD(5) * t70 - t44 * t65;
t86 = t45 ^ 2 - t48 ^ 2;
t51 = qJD(4) ^ 2;
t85 = t51 + t52;
t84 = qJD(2) * pkin(2);
t83 = qJD(2) * t18;
t82 = qJD(2) * t19;
t71 = pkin(4) * t80;
t68 = 0.2e1 * t72;
t28 = t95 * t48;
t15 = -pkin(7) * t80 + t48 * t24;
t13 = qJD(4) * pkin(4) + t15;
t67 = -pkin(4) * t41 - t13;
t66 = t85 * t49;
t39 = t45 * pkin(4) + qJ(3);
t64 = -t31 + t75;
t63 = -0.2e1 * t65;
t35 = pkin(4) * t78 + qJD(3);
t26 = -t44 * t45 + t90;
t21 = t39 * qJD(2) + t75;
t33 = t75 * t80;
t8 = t33 + (pkin(7) * qJD(2) - t24) * t79;
t9 = t24 * t78 + (-pkin(7) * t78 + t45 * t75) * qJD(2);
t59 = t21 * t18 - t44 * t9 + t47 * t8;
t57 = t26 * t41;
t56 = -t64 + t73;
t54 = t21 * t19 + (t96 * t14 - t8) * t44;
t29 = (qJD(3) + t74) * qJD(2);
t53 = t62 * qJD(2) - t50 * t51 + t29;
t30 = t62 - t84;
t27 = t95 * t45;
t23 = qJD(4) * t28;
t22 = t95 * t79;
t16 = (t35 + t74) * qJD(2);
t7 = t58 * t41;
t6 = qJD(2) * t61 + t87;
t3 = t18 ^ 2 - t19 ^ 2;
t2 = -t87 + (-t18 - t69) * t41;
t1 = t19 * t41 - t5;
t4 = [0, 0, -t89, -t88, t89, t88, t29 * t46 + (t92 + (t30 - t74) * t46) * qJD(2), 0, 0, 0, 0, 0, t48 * t46 * t68 + t45 * t66, t46 * t63 + t48 * t66, 0, 0, 0, 0, 0, (qJD(2) * t57 + t6) * t46 + ((t98 * t44 + t97 * t47) * t41 + t82) * t49, -0.2e1 * t46 * t5 + (-(-t98 * t47 + t55) * t41 - t83) * t49; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t29 * qJ(3) + t31 * qJD(3) + (-t92 + (-t30 - t84) * t46) * qJD(1), t48 * t63, t86 * t68, -t51 * t45, -t51 * t48, 0, t53 * t45 + t56 * t78, t53 * t48 - t56 * t79, t18 * t58 - t5 * t26, t18 * t11 + t19 * t58 + t5 * t25 - t26 * t6, -t7, -t94, 0, (t47 * t22 + t44 * t23 + (t27 * t47 + t28 * t44) * qJD(5)) * t41 + t35 * t19 + t39 * t6 + t16 * t25 + t21 * t11 + (-t49 * t19 - t46 * t57) * qJD(1), -(t44 * t22 - t47 * t23 + (t27 * t44 - t28 * t47) * qJD(5)) * t41 - t35 * t18 - t39 * t5 + t16 * t26 - t21 * t58 + (t49 * t18 + t46 * t58) * qJD(1); 0, 0, 0, 0, 0, -t52, t64 * qJD(2), 0, 0, 0, 0, 0, -t85 * t45, -t85 * t48, 0, 0, 0, 0, 0, -t7 - t82, t83 - t94; 0, 0, 0, 0, 0, 0, 0, t48 * t52 * t45, -t86 * t52, 0, 0, 0, -t31 * t80 + t33, -t64 * t81, -t93, t3, t1, t2, 0, -(-t44 * t15 - t91) * t41 - t19 * t71 + (t67 * t44 - t91) * qJD(5) + t59, t18 * t71 + (t67 * qJD(5) + t15 * t41 - t9) * t47 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t3, t1, t2, 0, t59 + t96 * (-t44 * t13 - t91), (-t96 * t13 - t9) * t47 + t54;];
tauc_reg = t4;
