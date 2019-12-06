% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRP1
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:39
% EndTime: 2019-12-05 15:28:41
% DurationCPUTime: 0.50s
% Computational Cost: add. (813->126), mult. (2132->156), div. (0->0), fcn. (1534->4), ass. (0->82)
t101 = cos(qJ(4));
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t63 = sin(qJ(4));
t44 = t101 * t61 + t63 * t62;
t105 = t44 * qJD(2);
t83 = t101 * t62;
t78 = qJD(2) * t83;
t51 = qJD(4) * t78;
t89 = qJD(2) * t61;
t84 = t63 * t89;
t27 = qJD(4) * t84 - t51;
t41 = t44 * qJD(4);
t96 = t63 * t61;
t72 = t83 - t96;
t108 = t105 * t41 + t27 * t72;
t28 = qJD(2) * t41;
t36 = -t78 + t84;
t80 = qJD(4) * t101;
t40 = qJD(4) * t96 - t62 * t80;
t92 = -t44 * t28 + t40 * t36;
t110 = t108 - t92;
t109 = t108 + t92;
t104 = t105 ^ 2;
t33 = t36 ^ 2;
t107 = -t33 - t104;
t106 = -t33 + t104;
t95 = pkin(6) + qJ(3);
t48 = t95 * t61;
t58 = t62 * qJD(1);
t31 = -qJD(2) * t48 + t58;
t85 = qJ(3) * qJD(2);
t46 = t61 * qJD(1) + t62 * t85;
t32 = t62 * qJD(2) * pkin(6) + t46;
t15 = t101 * t32 + t63 * t31;
t69 = t44 * qJD(3);
t67 = qJD(2) * t69;
t4 = t15 * qJD(4) + t67;
t49 = t95 * t62;
t73 = -t101 * t48 - t63 * t49;
t103 = t4 * t73;
t102 = t4 * t72;
t56 = -t62 * pkin(3) - pkin(2);
t47 = qJD(2) * t56 + qJD(3);
t10 = t36 * pkin(4) - qJ(5) * t105 + t47;
t100 = t10 * t105;
t99 = t105 * t36;
t97 = t63 * t32;
t11 = qJD(4) * qJ(5) + t15;
t94 = t11 - t15;
t91 = qJD(3) * t78 + t31 * t80;
t90 = t61 ^ 2 + t62 ^ 2;
t30 = qJD(4) * t41;
t12 = qJD(3) * t72 + qJD(4) * t73;
t88 = t12 * qJD(4);
t25 = t101 * t49 - t63 * t48;
t13 = qJD(4) * t25 + t69;
t87 = t13 * qJD(4);
t14 = t101 * t31 - t97;
t86 = qJD(5) - t14;
t82 = qJD(3) * t89;
t81 = t14 + t97;
t79 = t90 * qJD(2);
t77 = t28 * pkin(4) + t27 * qJ(5);
t75 = -t28 * t72 + t36 * t41;
t74 = (-t61 * t85 + t58) * t61 - t46 * t62;
t71 = t63 * t82 - t91;
t68 = t105 * t13 - t12 * t36 - t25 * t28 + t27 * t73 + t4 * t44;
t66 = 0.2e1 * t105 * qJD(4);
t65 = qJD(3) * t105;
t29 = t40 * qJD(4);
t19 = -pkin(4) * t72 - t44 * qJ(5) + t56;
t18 = pkin(4) * t105 + t36 * qJ(5);
t17 = t51 + (t36 - t84) * qJD(4);
t16 = -t51 + (t36 + t84) * qJD(4);
t8 = -qJD(4) * pkin(4) + t86;
t7 = t41 * pkin(4) + t40 * qJ(5) - t44 * qJD(5);
t6 = -qJD(5) * t105 + t77;
t3 = (-qJD(4) * t32 - t82) * t63 + t91;
t2 = -t105 * t40 - t27 * t44;
t1 = (qJD(5) - t97) * qJD(4) - t71;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, t109, -t14 * t41 - t15 * t40 + t3 * t44 - t102, 0, 0, 0, 0, 0, 0, -t30, t109, -t29, t1 * t44 - t11 * t40 + t8 * t41 - t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t79, (qJ(3) * t79 - t74) * qJD(3), t2, -t110, -t29, t75, -t30, 0, t56 * t28 + t47 * t41 - t87, -t56 * t27 - t47 * t40 - t88, t14 * t40 - t15 * t41 + t3 * t72 + t68, t15 * t12 - t14 * t13 + t3 * t25 - t103, t2, -t29, t110, 0, t30, t75, t10 * t41 + t19 * t28 + t7 * t36 - t6 * t72 - t87, t1 * t72 - t11 * t41 - t8 * t40 + t68, t10 * t40 - t105 * t7 + t19 * t27 - t6 * t44 + t88, t1 * t25 + t10 * t7 + t11 * t12 + t8 * t13 + t6 * t19 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 * qJD(2) ^ 2, t74 * qJD(2), 0, 0, 0, 0, 0, 0, t66, -t16, t107, t105 * t14 + t15 * t36, 0, 0, 0, 0, 0, 0, t66, t107, t16, t11 * t36 + (-qJD(5) - t8) * t105 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t106, t17, -t99, 0, 0, -t105 * t47 - t65, qJD(4) * t81 + t47 * t36 + t71, 0, 0, t99, t17, -t106, 0, 0, -t99, -t18 * t36 - t100 - t65, pkin(4) * t27 - t28 * qJ(5) + t94 * t105 + (t8 - t86) * t36, -t10 * t36 + t18 * t105 + (0.2e1 * qJD(5) - t81) * qJD(4) - t71, -t4 * pkin(4) + t1 * qJ(5) - t10 * t18 + t11 * t86 - t8 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t17, -qJD(4) ^ 2 - t104, -qJD(4) * t94 + t100 + t67;];
tauc_reg = t5;
