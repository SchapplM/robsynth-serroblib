% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP5
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
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:07
% DurationCPUTime: 0.49s
% Computational Cost: add. (811->131), mult. (2127->185), div. (0->0), fcn. (1332->4), ass. (0->89)
t86 = (qJD(1) * qJD(2));
t110 = -2 * t86;
t64 = sin(qJ(3));
t65 = sin(qJ(2));
t66 = cos(qJ(3));
t67 = cos(qJ(2));
t42 = t64 * t67 + t66 * t65;
t91 = qJD(1) * t42;
t61 = qJD(2) + qJD(3);
t109 = t91 ^ 2;
t108 = -pkin(6) - pkin(5);
t81 = qJD(2) * t108;
t44 = t65 * t81;
t46 = t67 * t81;
t48 = t108 * t65;
t49 = t108 * t67;
t74 = t66 * t48 + t64 * t49;
t5 = t74 * qJD(3) + t66 * t44 + t64 * t46;
t107 = t5 * t61;
t25 = t64 * t48 - t66 * t49;
t6 = t25 * qJD(3) + t64 * t44 - t66 * t46;
t106 = t6 * t61;
t45 = qJD(1) * t49;
t102 = t64 * t45;
t43 = qJD(1) * t48;
t92 = qJD(2) * pkin(2);
t38 = t43 + t92;
t18 = t66 * t38 + t102;
t105 = t18 * t61;
t100 = t66 * t45;
t19 = t64 * t38 - t100;
t104 = t19 * t61;
t99 = t66 * t67;
t82 = qJD(1) * t99;
t90 = qJD(1) * t65;
t83 = t64 * t90;
t34 = -t82 + t83;
t103 = t91 * t34;
t101 = t64 * t65;
t69 = qJD(1) ^ 2;
t98 = t67 * t69;
t68 = qJD(2) ^ 2;
t97 = t68 * t65;
t96 = t68 * t67;
t21 = t66 * t43 + t102;
t88 = qJD(3) * t66;
t95 = pkin(2) * t88 + qJD(4) - t21;
t80 = t67 * t86;
t94 = -qJD(3) * t82 - t66 * t80;
t93 = t65 ^ 2 - t67 ^ 2;
t89 = qJD(3) * t64;
t87 = qJD(4) - t18;
t85 = t65 * t92;
t84 = pkin(2) * t90;
t59 = -t67 * pkin(2) - pkin(1);
t76 = qJD(1) * t81;
t39 = t65 * t76;
t40 = t67 * t76;
t79 = -t38 * t88 - t66 * t39 - t64 * t40 - t45 * t89;
t4 = t38 * t89 + t64 * t39 - t66 * t40 - t45 * t88;
t78 = pkin(1) * t110;
t20 = t64 * t43 - t100;
t77 = pkin(2) * t89 - t20;
t47 = t59 * qJD(1);
t75 = t61 * t101;
t16 = pkin(3) * t91 + t34 * qJ(4);
t10 = t34 * pkin(3) - qJ(4) * t91 + t47;
t73 = -t10 * t91 - t4;
t72 = -t10 * t34 - t79;
t71 = -t47 * t91 - t4;
t70 = t47 * t34 + t79;
t23 = t61 * t42;
t60 = t61 * qJD(4);
t58 = -t66 * pkin(2) - pkin(3);
t56 = t64 * pkin(2) + qJ(4);
t41 = -t99 + t101;
t22 = -qJD(2) * t99 - t67 * t88 + t75;
t17 = t41 * pkin(3) - t42 * qJ(4) + t59;
t15 = t23 * qJD(1);
t14 = qJD(1) * t75 + t94;
t13 = t61 * qJ(4) + t19;
t12 = t16 + t84;
t11 = -t61 * pkin(3) + t87;
t9 = -t34 ^ 2 + t109;
t7 = -t94 + (t34 - t83) * t61;
t3 = t60 - t79;
t2 = t23 * pkin(3) + t22 * qJ(4) - t42 * qJD(4) + t85;
t1 = t15 * pkin(3) + t14 * qJ(4) + qJD(2) * t84 - qJD(4) * t91;
t8 = [0, 0, 0, 0.2e1 * t65 * t80, t93 * t110, t96, -t97, 0, -pkin(5) * t96 + t65 * t78, pkin(5) * t97 + t67 * t78, -t14 * t42 - t22 * t91, t14 * t41 - t42 * t15 + t22 * t34 - t23 * t91, -t22 * t61, -t23 * t61, 0, t59 * t15 + t47 * t23 - t106 + (qJD(1) * t41 + t34) * t85, -t59 * t14 - t47 * t22 + 0.2e1 * t91 * t85 - t107, t1 * t41 + t10 * t23 + t17 * t15 + t2 * t34 - t106, -t11 * t22 - t13 * t23 + t14 * t74 - t25 * t15 - t3 * t41 - t5 * t34 + t4 * t42 + t6 * t91, -t1 * t42 + t10 * t22 + t17 * t14 - t2 * t91 + t107, t1 * t17 + t10 * t2 + t11 * t6 + t13 * t5 + t3 * t25 - t4 * t74; 0, 0, 0, -t65 * t98, t93 * t69, 0, 0, 0, t69 * pkin(1) * t65, pkin(1) * t98, t103, t9, t7, 0, 0, t20 * t61 + (-t34 * t90 - t61 * t89) * pkin(2) + t71, t21 * t61 + (-t61 * t88 - t90 * t91) * pkin(2) + t70, -t12 * t34 - t77 * t61 + t73, -t58 * t14 - t56 * t15 + (t13 + t77) * t91 + (t11 - t95) * t34, t12 * t91 + t95 * t61 + t60 + t72, -t10 * t12 + t77 * t11 + t95 * t13 + t3 * t56 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t9, t7, 0, 0, t71 + t104, t70 + t105, -t16 * t34 + t104 + t73, pkin(3) * t14 - t15 * qJ(4) + (t13 - t19) * t91 + (t11 - t87) * t34, t16 * t91 - t105 + 0.2e1 * t60 + t72, -t4 * pkin(3) + t3 * qJ(4) - t10 * t16 - t11 * t19 + t87 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t7, -t61 ^ 2 - t109, -t13 * t61 - t73;];
tauc_reg = t8;
