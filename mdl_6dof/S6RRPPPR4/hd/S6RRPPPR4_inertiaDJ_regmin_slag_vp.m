% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:44
% EndTime: 2019-03-09 08:19:47
% DurationCPUTime: 1.03s
% Computational Cost: add. (817->166), mult. (1894->298), div. (0->0), fcn. (1471->6), ass. (0->97)
t121 = pkin(3) + pkin(7);
t79 = cos(qJ(6));
t104 = qJD(6) * t79;
t77 = sin(qJ(6));
t105 = qJD(6) * t77;
t74 = sin(pkin(9));
t75 = cos(pkin(9));
t34 = t74 * t104 - t75 * t105;
t78 = sin(qJ(2));
t107 = t78 * qJ(3);
t76 = -pkin(2) - qJ(4);
t80 = cos(qJ(2));
t110 = t76 * t80;
t120 = t107 - t110;
t48 = (t74 ^ 2 + t75 ^ 2) * qJD(4);
t108 = qJ(5) * t74;
t117 = pkin(4) + pkin(5);
t119 = -t117 * t75 - t108;
t118 = 2 * qJD(3);
t116 = pkin(8) + t76;
t58 = -t75 * qJD(5) + qJD(3);
t115 = t58 * t80;
t114 = t74 * t78;
t113 = t74 * t80;
t112 = t75 * t78;
t111 = t75 * t80;
t103 = t78 * qJD(2);
t93 = pkin(2) * t103 - t78 * qJD(3);
t27 = -t80 * qJD(4) + (-qJ(3) * t80 + qJ(4) * t78) * qJD(2) + t93;
t68 = t80 * qJD(2);
t64 = pkin(7) * t68;
t44 = pkin(3) * t68 + t64;
t13 = t75 * t27 + t74 * t44;
t38 = -pkin(1) - t120;
t52 = t121 * t78;
t23 = t75 * t38 + t74 * t52;
t53 = t121 * t80;
t106 = qJD(4) * t78;
t102 = t80 * qJD(3);
t101 = qJ(3) * qJD(3);
t100 = qJD(2) * qJ(3);
t99 = -0.2e1 * pkin(1) * qJD(2);
t20 = t78 * qJ(5) + t23;
t96 = t74 * t68;
t95 = t75 * t103;
t56 = t75 * t68;
t12 = -t74 * t27 + t75 * t44;
t22 = -t74 * t38 + t75 * t52;
t7 = qJ(5) * t68 + t78 * qJD(5) + t13;
t94 = t75 * qJ(5) - qJ(3);
t92 = t76 * t48;
t63 = pkin(7) * t103;
t91 = -qJD(5) * t113 + t63;
t8 = -pkin(4) * t68 - t12;
t3 = t7 * t74 - t8 * t75;
t90 = -t80 * pkin(2) - t107;
t89 = pkin(4) * t75 + t108;
t11 = pkin(8) * t113 - t117 * t78 - t22;
t15 = pkin(8) * t111 + t20;
t88 = t79 * t11 - t77 * t15;
t87 = t77 * t11 + t79 * t15;
t4 = t12 * t75 + t13 * t74;
t45 = t116 * t74;
t46 = t116 * t75;
t86 = t79 * t45 - t77 * t46;
t85 = t77 * t45 + t79 * t46;
t84 = t79 * t74 - t77 * t75;
t83 = t77 * t74 + t79 * t75;
t82 = -t34 * t78 - t68 * t83;
t35 = t83 * qJD(6);
t81 = t90 * qJD(2) + t102;
t55 = t74 * t103;
t54 = 0.2e1 * t78 * t68;
t50 = -pkin(1) + t90;
t49 = t76 * t56;
t47 = t74 * pkin(4) - t94;
t43 = -pkin(3) * t103 - t63;
t39 = 0.2e1 * t48;
t36 = -t117 * t74 + t94;
t33 = -t80 * t100 + t93;
t30 = t84 * t80;
t29 = t83 * t80;
t28 = t89 * t80 + t53;
t26 = t119 * t80 - t53;
t21 = -t78 * pkin(4) - t22;
t19 = (-pkin(3) - t89) * t103 - t91;
t18 = -t35 * t78 + t68 * t84;
t17 = t84 * t103 + t80 * t35;
t16 = t83 * t103 - t34 * t80;
t14 = (pkin(3) - t119) * t103 + t91;
t10 = t83 * qJD(4) - t86 * qJD(6);
t9 = t84 * qJD(4) + t85 * qJD(6);
t6 = -pkin(8) * t95 + t7;
t5 = (-pkin(8) * t114 - t117 * t80) * qJD(2) - t12;
t2 = -t87 * qJD(6) + t79 * t5 - t77 * t6;
t1 = -t88 * qJD(6) - t77 * t5 - t79 * t6;
t24 = [0, 0, 0, t54, 0.2e1 * (-t78 ^ 2 + t80 ^ 2) * qJD(2), 0, 0, 0, t78 * t99, t80 * t99, 0, -0.2e1 * t50 * t103 + 0.2e1 * t33 * t80, -0.2e1 * t33 * t78 - 0.2e1 * t50 * t68, 0.2e1 * t50 * t33, 0.2e1 * t43 * t111 + 0.2e1 * t12 * t78 + 0.2e1 * (-t53 * t112 + t22 * t80) * qJD(2), -0.2e1 * t43 * t113 - 0.2e1 * t13 * t78 + 0.2e1 * (t53 * t114 - t23 * t80) * qJD(2), 0.2e1 * (t12 * t74 - t13 * t75) * t80 + 0.2e1 * (-t22 * t74 + t23 * t75) * t103, 0.2e1 * t22 * t12 + 0.2e1 * t23 * t13 + 0.2e1 * t53 * t43, 0.2e1 * t19 * t111 - 0.2e1 * t8 * t78 + 0.2e1 * (-t28 * t112 - t21 * t80) * qJD(2), 0.2e1 * (-t7 * t75 - t74 * t8) * t80 + 0.2e1 * (t20 * t75 + t21 * t74) * t103, 0.2e1 * t19 * t113 + 0.2e1 * t7 * t78 + 0.2e1 * (-t28 * t114 + t20 * t80) * qJD(2), 0.2e1 * t28 * t19 + 0.2e1 * t20 * t7 + 0.2e1 * t21 * t8, -0.2e1 * t30 * t17, 0.2e1 * t30 * t16 + 0.2e1 * t17 * t29, -0.2e1 * t17 * t78 + 0.2e1 * t30 * t68, 0.2e1 * t16 * t78 - 0.2e1 * t29 * t68, t54, -0.2e1 * t14 * t29 + 0.2e1 * t26 * t16 - 0.2e1 * t2 * t78 - 0.2e1 * t88 * t68, -0.2e1 * t1 * t78 - 0.2e1 * t14 * t30 + 0.2e1 * t26 * t17 + 0.2e1 * t87 * t68; 0, 0, 0, 0, 0, t68, -t103, 0, -t64, t63, t81, t64, -t63, t81 * pkin(7), t43 * t74 + t49 + (t102 + (-qJD(4) - t100) * t78) * t75, t43 * t75 + (t120 * qJD(2) - t102 + t106) * t74, -t4, t43 * qJ(3) + t53 * qJD(3) + t4 * t76 + (-t22 * t75 - t23 * t74) * qJD(4), t19 * t74 + t49 + (t115 + (-qJD(2) * t47 - qJD(4)) * t78) * t75, -t3, -t19 * t75 + (-t106 + t115 + (-t47 * t78 + t110) * qJD(2)) * t74, t19 * t47 + t28 * t58 + t3 * t76 + (-t20 * t74 + t21 * t75) * qJD(4), t17 * t83 - t30 * t34, -t16 * t83 + t17 * t84 + t34 * t29 + t30 * t35, t82, -t18, 0, -t10 * t78 - t14 * t84 + t36 * t16 + t26 * t35 + t58 * t29 + t85 * t68, t14 * t83 + t36 * t17 + t26 * t34 + t58 * t30 + t86 * t68 - t9 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0.2e1 * t101, t74 * t118, t75 * t118, t39, -0.2e1 * t92 + 0.2e1 * t101, 0.2e1 * t58 * t74, t39, -0.2e1 * t58 * t75, 0.2e1 * t47 * t58 - 0.2e1 * t92, 0.2e1 * t83 * t34, 0.2e1 * t34 * t84 - 0.2e1 * t35 * t83, 0, 0, 0, 0.2e1 * t36 * t35 + 0.2e1 * t58 * t84, 0.2e1 * t36 * t34 - 0.2e1 * t58 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, t64, t56, -t96, 0, t4, t56, 0, t96, t3, 0, 0, 0, 0, 0, -t82, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, -t48, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t55, 0, t43, -t95, 0, -t55, t19, 0, 0, 0, 0, 0, -t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, t58, 0, 0, 0, 0, 0, -t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t55, 0, t8, 0, 0, 0, 0, 0, t78 * t105 - t79 * t68, t78 * t104 + t77 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t68, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t35, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t24;
