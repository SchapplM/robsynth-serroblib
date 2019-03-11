% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:14
% EndTime: 2019-03-09 01:30:16
% DurationCPUTime: 0.70s
% Computational Cost: add. (571->146), mult. (1182->217), div. (0->0), fcn. (673->6), ass. (0->90)
t43 = sin(qJ(5));
t79 = t43 * qJD(1);
t31 = qJD(6) + t79;
t108 = qJD(6) - t31;
t32 = sin(pkin(9)) * pkin(1) + qJ(3);
t107 = qJD(1) * t32;
t45 = cos(qJ(5));
t42 = sin(qJ(6));
t81 = qJD(6) * t42;
t70 = t45 * t81;
t44 = cos(qJ(6));
t78 = t44 * qJD(5);
t39 = t45 ^ 2;
t86 = qJD(1) * t39;
t106 = -(-t31 * t43 + t86) * t78 + t31 * t70;
t85 = qJD(1) * t45;
t68 = t44 * t85;
t23 = t42 * qJD(5) + t68;
t76 = qJD(1) * qJD(5);
t65 = t43 * t76;
t10 = qJD(6) * t23 - t42 * t65;
t25 = qJD(4) + t107;
t17 = -qJD(1) * pkin(7) + t25;
t12 = t45 * qJD(2) + t43 * t17;
t75 = qJD(3) * qJD(1);
t4 = t12 * qJD(5) - t45 * t75;
t105 = t4 * t42;
t104 = t4 * t44;
t49 = -t43 * t78 - t70;
t9 = qJD(1) * t49 + qJD(6) * t78;
t102 = t9 * t42;
t82 = qJD(5) * t45;
t101 = t23 * t82 + t9 * t43;
t21 = t42 * t85 - t78;
t100 = t21 * t31;
t99 = t23 * t31;
t98 = t42 * t31;
t97 = t42 * t43;
t96 = t43 * t10;
t95 = t44 * t31;
t94 = t44 * t45;
t93 = t45 * t21;
t92 = t45 * t23;
t46 = qJD(5) ^ 2;
t91 = t46 * t43;
t90 = t46 * t45;
t89 = t43 ^ 2 - t39;
t47 = qJD(1) ^ 2;
t88 = -t46 - t47;
t30 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t87 = qJD(1) * t30;
t29 = -pkin(7) + t32;
t84 = qJD(5) * t29;
t83 = qJD(5) * t43;
t80 = qJD(6) * t44;
t18 = -qJD(3) + t87;
t77 = qJD(3) - t18;
t74 = t31 * t97;
t73 = t43 * t95;
t72 = t42 * t86;
t71 = t31 * t81;
t69 = t31 * t80;
t67 = 0.2e1 * qJD(4) * qJD(1);
t7 = qJD(5) * pkin(8) + t12;
t66 = t29 * t31 + t7;
t64 = t45 * t76;
t63 = t31 + t79;
t62 = t77 * qJD(1);
t60 = -0.2e1 * t64;
t59 = qJD(6) * t43 + qJD(1);
t58 = pkin(5) * t45 + pkin(8) * t43;
t57 = t42 * t64 + t69;
t16 = t43 * pkin(5) - t45 * pkin(8) + t30;
t8 = qJD(1) * t16 - qJD(3);
t2 = t42 * t8 + t44 * t7;
t56 = t42 * t7 - t44 * t8;
t55 = qJD(3) + t18 + t87;
t53 = t43 * qJD(2) - t45 * t17;
t52 = -t29 * t46 + t67;
t51 = qJD(5) * t74 - t45 * t69;
t6 = -qJD(5) * pkin(5) + t53;
t50 = -pkin(8) * t82 + t43 * t6;
t20 = t58 * qJD(5) + qJD(4);
t3 = -qJD(5) * t53 + t43 * t75;
t48 = -qJD(3) * t31 - t6 * qJD(5) - qJD(6) * t8 - t3;
t36 = 0.2e1 * t75;
t24 = t58 * qJD(1);
t14 = t20 * qJD(1);
t13 = t44 * t14;
t1 = [0, 0, 0, 0, 0, t36, 0.2e1 * t107 * qJD(3), t36, t67, t25 * qJD(3) + t18 * qJD(4) + (qJD(3) * t32 + qJD(4) * t30) * qJD(1), t43 * t60, 0.2e1 * t89 * t76, -t91, -t90, 0, t43 * t52 + t55 * t82, t45 * t52 - t55 * t83, t23 * t49 + t9 * t94 (t21 * t44 + t23 * t42) * t83 + (-t10 * t44 - t102 + (t21 * t42 - t23 * t44) * qJD(6)) * t45, t101 - t106, -t96 + (-t72 - t93) * qJD(5) + t51, t63 * t82 (-t16 * t81 + t44 * t20) * t31 + (t21 * t84 + t48 * t42 - t66 * t80 + t13) * t43 + (t6 * t80 - qJD(3) * t21 - t29 * t10 + t105 + (-t29 * t98 + (t44 * t16 - t29 * t97) * qJD(1) - t56) * qJD(5)) * t45 -(t16 * t80 + t42 * t20) * t31 + (t23 * t84 + (t66 * qJD(6) - t14) * t42 + t48 * t44) * t43 + (-t6 * t81 - qJD(3) * t23 - t29 * t9 + t104 + (-t29 * t95 - (t44 * t43 * t29 + t42 * t16) * qJD(1) - t2) * qJD(5)) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t91, 0, 0, 0, 0, 0, t96 + (-t72 + t93) * qJD(5) + t51, t101 + t106; 0, 0, 0, 0, 0, -t47, -t107 * qJD(1), -t47, 0 (-qJD(4) - t25) * qJD(1), 0, 0, 0, 0, 0, t60, 0.2e1 * t65, 0, 0, 0, 0, 0, t71 + (t74 + (t21 - t78) * t45) * qJD(1) (t73 + t92) * qJD(1) + t57; 0, 0, 0, 0, 0, 0, 0, 0, -t47, t62, 0, 0, 0, 0, 0, t88 * t43, t88 * t45, 0, 0, 0, 0, 0, -t45 * t10 - t59 * t95 + (-t63 * t45 * t42 + t21 * t43) * qJD(5), -t45 * t9 + t59 * t98 + (-t31 * t94 + (t23 - t68) * t43) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t47 * t43, -t89 * t47, 0, 0, 0, t45 * t62, -t77 * t79, t23 * t95 + t102 (t9 - t100) * t44 + (-t10 - t99) * t42 (t73 - t92) * qJD(1) + t57, -t71 + (-t74 + (t21 + t78) * t45) * qJD(1), -t31 * t85, -pkin(5) * t10 - t104 - (t44 * t24 + t42 * t53) * t31 - t12 * t21 + (-pkin(8) * t95 + t6 * t42) * qJD(6) + (t42 * t50 + t45 * t56) * qJD(1), -pkin(5) * t9 + t105 + (t42 * t24 - t44 * t53) * t31 - t12 * t23 + (pkin(8) * t98 + t6 * t44) * qJD(6) + (t2 * t45 + t44 * t50) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t21, -t21 ^ 2 + t23 ^ 2, t9 + t100, -t10 + t99, t64, -t108 * t2 - t6 * t23 - t42 * t3 + t13, t108 * t56 - t42 * t14 + t6 * t21 - t44 * t3;];
tauc_reg  = t1;
