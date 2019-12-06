% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:41
% EndTime: 2019-12-05 15:38:44
% DurationCPUTime: 0.51s
% Computational Cost: add. (731->134), mult. (1910->177), div. (0->0), fcn. (1351->6), ass. (0->84)
t59 = cos(pkin(8));
t62 = cos(qJ(4));
t95 = t62 * t59;
t58 = sin(pkin(8));
t60 = sin(qJ(4));
t98 = t60 * t58;
t43 = -t95 + t98;
t44 = t62 * t58 + t60 * t59;
t63 = cos(qJ(2));
t85 = t63 * qJD(1);
t74 = qJD(3) - t85;
t61 = sin(qJ(2));
t86 = t61 * qJD(1);
t50 = qJD(2) * qJ(3) + t86;
t91 = t58 ^ 2 + t59 ^ 2;
t106 = t91 * t50;
t41 = t44 * qJD(4);
t69 = t43 * t63;
t92 = pkin(6) + qJ(3);
t46 = t92 * t58;
t47 = t92 * t59;
t73 = -t62 * t46 - t60 * t47;
t101 = -qJD(1) * t69 + t43 * qJD(3) - t73 * qJD(4);
t105 = t101 * qJD(4);
t67 = qJD(2) * t44;
t104 = t67 ^ 2;
t82 = qJD(2) * t95;
t83 = qJD(2) * t98;
t36 = -t82 + t83;
t54 = -t59 * pkin(3) - pkin(2);
t42 = t54 * qJD(2) + t74;
t5 = t36 * pkin(4) - qJ(5) * t67 + t42;
t103 = t5 * t67;
t19 = -t60 * t46 + t62 * t47;
t102 = -t19 * qJD(4) - t74 * t44;
t100 = t67 * t36;
t77 = pkin(6) * qJD(2) + t50;
t27 = t77 * t59;
t99 = t60 * t27;
t64 = qJD(2) ^ 2;
t94 = t64 * t61;
t93 = t64 * t63;
t90 = qJD(2) * pkin(2);
t89 = qJD(2) * t61;
t88 = qJD(4) * t60;
t87 = qJD(4) * t62;
t26 = t77 * t58;
t10 = -t62 * t26 - t99;
t84 = qJD(5) - t10;
t81 = t91 * t63;
t45 = (qJD(3) + t85) * qJD(2);
t80 = t91 * t45;
t79 = t10 + t99;
t78 = t102 * qJD(4);
t40 = t58 * t88 - t59 * t87;
t76 = -t41 * pkin(4) - t40 * qJ(5) + t44 * qJD(5) + t86;
t75 = t91 * qJD(3);
t2 = -t26 * t88 + t27 * t87 + t44 * t45;
t11 = -t60 * t26 + t62 * t27;
t72 = t26 * t87 + t43 * t45;
t49 = qJD(4) * t82;
t28 = qJD(4) * t83 - t49;
t29 = qJD(2) * t41;
t55 = qJD(2) * t86;
t71 = t29 * pkin(4) + t28 * qJ(5) + t55;
t31 = t43 * t61;
t13 = -qJD(4) * t31 + t63 * t67;
t70 = -t13 * qJD(4) - t63 * t29 + t36 * t89;
t68 = t11 * qJD(4) - t2;
t12 = -qJD(2) * t69 - t61 * t41;
t66 = t12 * qJD(4) - t63 * t28 - t67 * t89;
t65 = 0.2e1 * t67 * qJD(4);
t48 = t74 - t90;
t32 = t36 ^ 2;
t30 = t44 * t61;
t17 = t43 * pkin(4) - t44 * qJ(5) + t54;
t16 = pkin(4) * t67 + t36 * qJ(5);
t15 = t49 + (t36 - t83) * qJD(4);
t14 = -t49 + (t36 + t83) * qJD(4);
t7 = qJD(4) * qJ(5) + t11;
t6 = -qJD(4) * pkin(4) + t84;
t3 = -qJD(5) * t67 + t71;
t1 = (qJD(5) - t99) * qJD(4) - t72;
t4 = [0, 0, -t94, -t93, -t59 * t94, t58 * t94, t91 * t93, t61 * t80 + (t48 * t61 + (-t86 + t106) * t63) * qJD(2), 0, 0, 0, 0, 0, t70, -t66, t70, -t12 * t36 + t13 * t67 - t30 * t28 + t31 * t29, t66, -t1 * t31 + t7 * t12 + t6 * t13 + t2 * t30 - t3 * t63 + t5 * t89; 0, 0, 0, 0, 0, 0, t80 + (-qJD(1) * t81 + t75) * qJD(2), t50 * t75 + qJ(3) * t80 + ((-t48 - t90) * t61 - t50 * t81) * qJD(1), -t28 * t44 - t40 * t67, t28 * t43 - t44 * t29 + t40 * t36 - t41 * t67, -t40 * qJD(4), -t41 * qJD(4), 0, t54 * t29 + t42 * t41 + t78 + (qJD(2) * t43 - t36) * t86, -t54 * t28 - t42 * t40 + t105, t17 * t29 + t3 * t43 - t76 * t36 + t5 * t41 + t78, -t1 * t43 + t101 * t36 - t102 * t67 - t19 * t29 + t2 * t44 + t28 * t73 - t6 * t40 - t7 * t41, t17 * t28 - t3 * t44 + t5 * t40 + t67 * t76 - t105, t1 * t19 - t101 * t7 - t102 * t6 + t3 * t17 - t2 * t73 - t76 * t5; 0, 0, 0, 0, 0, 0, -t91 * t64, -qJD(2) * t106 + t55, 0, 0, 0, 0, 0, t65, -t14, t65, -t32 - t104, t14, t7 * t36 + (-qJD(5) - t6) * t67 + t71; 0, 0, 0, 0, 0, 0, 0, 0, t100, -t32 + t104, t15, 0, 0, -t42 * t67 + t68, t79 * qJD(4) + t42 * t36 + t72, -t16 * t36 - t103 + t68, pkin(4) * t28 - t29 * qJ(5) + (-t11 + t7) * t67 + (t6 - t84) * t36, t16 * t67 - t5 * t36 + (0.2e1 * qJD(5) - t79) * qJD(4) - t72, -t2 * pkin(4) + t1 * qJ(5) - t6 * t11 - t5 * t16 + t84 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t15, -qJD(4) ^ 2 - t104, -t7 * qJD(4) + t103 + t2;];
tauc_reg = t4;
