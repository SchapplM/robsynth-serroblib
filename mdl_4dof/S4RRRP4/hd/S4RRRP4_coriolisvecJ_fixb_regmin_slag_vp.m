% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP4
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
% tauc_reg [4x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:43
% DurationCPUTime: 0.51s
% Computational Cost: add. (630->117), mult. (1723->173), div. (0->0), fcn. (1112->4), ass. (0->89)
t108 = pkin(5) + pkin(6);
t63 = cos(qJ(2));
t47 = t108 * t63;
t43 = qJD(1) * t47;
t60 = sin(qJ(3));
t29 = t60 * t43;
t61 = sin(qJ(2));
t46 = t108 * t61;
t41 = qJD(1) * t46;
t95 = qJD(2) * pkin(2);
t35 = -t41 + t95;
t62 = cos(qJ(3));
t79 = t62 * t35 - t29;
t40 = t60 * t63 + t62 * t61;
t92 = qJD(1) * t40;
t93 = t92 * qJ(4);
t112 = t93 - t79;
t88 = qJD(1) * qJD(2);
t111 = -0.2e1 * t88;
t81 = qJD(2) * t108;
t73 = qJD(1) * t81;
t36 = t61 * t73;
t110 = (qJD(3) * t35 - t36) * t62;
t57 = qJD(2) + qJD(3);
t109 = t92 ^ 2;
t5 = t57 * pkin(3) - t112;
t107 = t5 + t112;
t102 = t62 * t63;
t84 = qJD(1) * t102;
t91 = qJD(1) * t61;
t85 = t60 * t91;
t25 = -t84 + t85;
t56 = -t63 * pkin(2) - pkin(1);
t45 = t56 * qJD(1);
t19 = t25 * pkin(3) + qJD(4) + t45;
t106 = t19 * t92;
t105 = t92 * t25;
t104 = t45 * t92;
t103 = t60 * t61;
t33 = t62 * t43;
t65 = qJD(1) ^ 2;
t101 = t63 * t65;
t64 = qJD(2) ^ 2;
t100 = t64 * t61;
t99 = t64 * t63;
t98 = -t62 * t41 - t29;
t80 = t63 * t88;
t97 = -qJD(3) * t84 - t62 * t80;
t96 = t61 ^ 2 - t63 ^ 2;
t94 = t25 * qJ(4);
t90 = qJD(3) * t60;
t89 = qJD(3) * t62;
t87 = t61 * t95;
t86 = pkin(2) * t91;
t83 = -pkin(2) * t57 - t35;
t18 = t57 * t40;
t16 = t18 * qJD(1);
t82 = t16 * pkin(3) + qJD(2) * t86;
t37 = t63 * t73;
t78 = t60 * t36 - t62 * t37;
t77 = -t60 * t37 - t43 * t90;
t76 = t60 * t41 - t33;
t74 = pkin(1) * t111;
t72 = t57 * t103;
t71 = -t60 * t35 - t33;
t70 = t60 * t46 - t62 * t47;
t69 = t45 * t25 - t77;
t42 = t61 * t81;
t44 = t63 * t81;
t68 = -t62 * t42 - t60 * t44 - t46 * t89 - t47 * t90;
t67 = t71 * qJD(3) + t78;
t66 = t70 * qJD(3) + t60 * t42 - t62 * t44;
t55 = t62 * pkin(2) + pkin(3);
t39 = -t102 + t103;
t24 = t25 ^ 2;
t17 = -qJD(2) * t102 - t63 * t89 + t72;
t15 = qJD(1) * t72 + t97;
t14 = -t39 * qJ(4) - t70;
t13 = -t40 * qJ(4) - t62 * t46 - t60 * t47;
t12 = -t24 + t109;
t11 = -t93 + t98;
t10 = t76 + t94;
t9 = -t71 - t94;
t6 = -t97 + (t25 - t85) * t57;
t4 = t17 * qJ(4) - t40 * qJD(4) + t66;
t3 = -t18 * qJ(4) - t39 * qJD(4) + t68;
t2 = t15 * qJ(4) - qJD(4) * t92 + t67;
t1 = -t16 * qJ(4) - t25 * qJD(4) + t110 + t77;
t7 = [0, 0, 0, 0.2e1 * t61 * t80, t96 * t111, t99, -t100, 0, -pkin(5) * t99 + t61 * t74, pkin(5) * t100 + t63 * t74, -t15 * t40 - t17 * t92, t15 * t39 - t40 * t16 + t17 * t25 - t18 * t92, -t17 * t57, -t18 * t57, 0, t56 * t16 + t45 * t18 + t66 * t57 + (qJD(1) * t39 + t25) * t87, -t56 * t15 - t45 * t17 - t68 * t57 + 0.2e1 * t92 * t87, -t1 * t39 + t13 * t15 - t14 * t16 + t5 * t17 - t9 * t18 - t2 * t40 - t3 * t25 - t4 * t92, t1 * t14 + t9 * t3 + t2 * t13 + t5 * t4 + t82 * (t39 * pkin(3) + t56) + t19 * (t18 * pkin(3) + t87); 0, 0, 0, -t61 * t101, t96 * t65, 0, 0, 0, t65 * pkin(1) * t61, pkin(1) * t101, t105, t12, t6, 0, 0, -t25 * t86 - t104 - t76 * t57 + (t83 * t60 - t33) * qJD(3) + t78, -t92 * t86 + t98 * t57 + (t83 * qJD(3) + t36) * t62 + t69, t55 * t15 + (t10 + t9) * t92 + (t11 - t5) * t25 + (-t16 * t60 + (-t25 * t62 + t60 * t92) * qJD(3)) * pkin(2), -pkin(3) * t106 - t5 * t10 - t9 * t11 + t2 * t55 + (-t19 * t91 + t1 * t60 + (-t5 * t60 + t62 * t9) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t12, t6, 0, 0, -t71 * t57 - t104 + t67, t79 * t57 - t110 + t69, pkin(3) * t15 - t107 * t25, t107 * t9 + (t2 - t106) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 - t109, t9 * t25 + t5 * t92 + t82;];
tauc_reg = t7;
