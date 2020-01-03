% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:52
% DurationCPUTime: 0.50s
% Computational Cost: add. (911->134), mult. (1881->189), div. (0->0), fcn. (1099->4), ass. (0->86)
t74 = -pkin(1) - pkin(6);
t52 = t74 * qJD(1) + qJD(2);
t120 = -qJ(4) * qJD(1) + t52;
t70 = sin(pkin(7));
t71 = cos(pkin(7));
t72 = sin(qJ(3));
t73 = cos(qJ(3));
t85 = t70 * t73 + t71 * t72;
t118 = t85 * qJD(1);
t33 = qJD(3) * t118;
t101 = qJD(1) * t73;
t102 = qJD(1) * t72;
t88 = t70 * t102;
t43 = t71 * t101 - t88;
t38 = t43 ^ 2;
t119 = -t118 ^ 2 - t38;
t103 = qJ(4) - t74;
t86 = t103 * t73;
t97 = t72 * qJD(4);
t34 = -qJD(3) * t86 - t97;
t100 = qJD(3) * t72;
t96 = t73 * qJD(4);
t80 = t103 * t100 - t96;
t13 = t70 * t34 - t71 * t80;
t14 = t71 * t34 + t70 * t80;
t50 = t103 * t72;
t24 = -t70 * t50 + t71 * t86;
t25 = -t71 * t50 - t70 * t86;
t90 = qJD(1) * qJD(3);
t87 = t73 * t90;
t51 = t71 * t87;
t32 = qJD(3) * t88 - t51;
t117 = -t118 * t14 + t13 * t43 - t24 * t33 + t25 * t32;
t66 = qJD(1) * qJD(2);
t116 = 0.2e1 * t66;
t91 = qJ(4) * qJD(3);
t99 = qJD(3) * t73;
t26 = t52 * t99 + (-t73 * t91 - t97) * qJD(1);
t77 = -t52 * t100 + (t72 * t91 - t96) * qJD(1);
t3 = t70 * t26 - t71 * t77;
t115 = t3 * t24;
t46 = -t70 * t72 + t71 * t73;
t114 = t3 * t46;
t35 = t120 * t72;
t111 = t70 * t35;
t29 = t71 * t35;
t75 = qJD(3) ^ 2;
t110 = t75 * t72;
t109 = t75 * t73;
t4 = t71 * t26 + t70 * t77;
t36 = t120 * t73;
t31 = qJD(3) * pkin(3) + t36;
t11 = t70 * t31 + t29;
t108 = pkin(3) * t87 + t66;
t107 = t72 ^ 2 - t73 ^ 2;
t76 = qJD(1) ^ 2;
t106 = -t75 - t76;
t105 = t76 * qJ(2);
t104 = t72 * pkin(3) + qJ(2);
t49 = pkin(3) * t102 + qJD(1) * qJ(2) + qJD(4);
t98 = t49 * qJD(1);
t95 = pkin(3) * t99 + qJD(2);
t16 = t71 * t36 - t111;
t94 = qJD(5) - t16;
t93 = qJ(2) * qJD(3);
t89 = 0.2e1 * qJD(1);
t10 = t71 * t31 - t111;
t12 = pkin(4) * t118 - t43 * qJ(5) + t49;
t84 = t12 * t43 + t3;
t82 = -t32 * pkin(4) + t33 * qJ(5) + t108;
t41 = t70 * t100 - t71 * t99;
t42 = -t71 * t100 - t70 * t99;
t81 = t118 * t41 + t32 * t85 + t46 * t33 - t42 * t43;
t2 = qJD(3) * qJD(5) + t4;
t7 = -qJD(3) * pkin(4) + qJD(5) - t10;
t8 = qJD(3) * qJ(5) + t11;
t79 = t2 * t85 - t8 * t41 - t7 * t42 - t114;
t78 = t10 * t42 - t11 * t41 + t4 * t85 - t114;
t62 = -t71 * pkin(3) - pkin(4);
t58 = t70 * pkin(3) + qJ(5);
t18 = pkin(4) * t85 - t46 * qJ(5) + t104;
t17 = pkin(3) * t101 + t43 * pkin(4) + qJ(5) * t118;
t15 = t70 * t36 + t29;
t6 = -t41 * pkin(4) - t42 * qJ(5) - t46 * qJD(5) + t95;
t1 = -t43 * qJD(5) + t82;
t5 = [0, 0, 0, 0, t116, qJ(2) * t116, -0.2e1 * t72 * t87, 0.2e1 * t107 * t90, -t110, -t109, 0, -t74 * t110 + (qJD(2) * t72 + t73 * t93) * t89, -t74 * t109 + (qJD(2) * t73 - t72 * t93) * t89, -t78 + t117, -t10 * t13 + t108 * t104 + t11 * t14 + t4 * t25 + t49 * t95 + t115, -t13 * qJD(3) + t1 * t85 + t118 * t6 - t12 * t41 - t18 * t32, -t79 + t117, t14 * qJD(3) - t1 * t46 - t12 * t42 + t18 * t33 - t6 * t43, t1 * t18 + t12 * t6 + t7 * t13 + t8 * t14 + t2 * t25 + t115; 0, 0, 0, 0, -t76, -t105, 0, 0, 0, 0, 0, t106 * t72, t106 * t73, t81, t78 - t98, -qJD(1) * t118 + t42 * qJD(3), t81, qJD(1) * t43 - t41 * qJD(3), -t12 * qJD(1) + t79; 0, 0, 0, 0, 0, 0, t73 * t76 * t72, -t107 * t76, 0, 0, 0, -t73 * t105, t72 * t105, (t11 - t15) * t43 + (-t10 + t16) * t118 + (t32 * t70 + t33 * t71) * pkin(3), t10 * t15 - t11 * t16 + (-t3 * t71 + t4 * t70 - t73 * t98) * pkin(3), t15 * qJD(3) - t118 * t17 - t84, t58 * t32 - t62 * t33 + (-t15 + t8) * t43 + (t7 - t94) * t118, -t12 * t118 + t17 * t43 + (0.2e1 * qJD(5) - t16) * qJD(3) + t4, -t12 * t17 - t7 * t15 + t2 * t58 + t3 * t62 + t94 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t10 * t43 + t11 * t118 + t108, t51 + (t43 - t88) * qJD(3), t119, 0.2e1 * t33, t8 * t118 + (-qJD(5) - t7) * t43 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t118, 0, -t38 - t75, -t8 * qJD(3) + t84;];
tauc_reg = t5;
