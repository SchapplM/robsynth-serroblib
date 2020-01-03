% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:31
% DurationCPUTime: 0.45s
% Computational Cost: add. (934->109), mult. (1456->143), div. (0->0), fcn. (636->4), ass. (0->85)
t71 = qJD(1) - qJD(3);
t105 = t71 ^ 2;
t46 = -pkin(1) - pkin(2);
t32 = t46 * qJD(1) + qJD(2);
t72 = (qJD(1) * qJD(2));
t104 = qJD(3) * t32 + t72;
t103 = qJD(4) * t71;
t45 = cos(qJ(3));
t102 = t45 * t105;
t101 = 0.2e1 * t45 * t103;
t43 = sin(qJ(3));
t74 = qJD(1) * qJ(2);
t63 = qJD(3) * t74;
t14 = t104 * t43 + t45 * t63;
t54 = t45 * qJ(2) + t43 * t46;
t21 = t43 * qJD(2) + t54 * qJD(3);
t100 = t21 * t71 + t14;
t99 = t43 * qJ(2) - t45 * t46;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t55 = pkin(4) * t42 - qJ(5) * t44;
t25 = t55 * qJD(4) - t42 * qJD(5);
t13 = t104 * t45 - t43 * t63;
t11 = t44 * t13;
t23 = t43 * t32 + t45 * t74;
t17 = -pkin(7) * t71 + t23;
t87 = t42 * t17;
t1 = t11 + (qJD(5) - t87) * qJD(4);
t10 = t42 * t13;
t77 = qJD(4) * t44;
t2 = t17 * t77 + t10;
t62 = -qJD(4) * pkin(4) + qJD(5);
t8 = t62 + t87;
t73 = qJD(4) * qJ(5);
t9 = t44 * t17 + t73;
t57 = t42 * t9 - t44 * t8;
t49 = -t57 * qJD(4) + t1 * t44 + t2 * t42;
t47 = qJD(4) ^ 2;
t98 = pkin(7) * t47;
t97 = t71 * pkin(3);
t20 = t45 * qJD(2) - t99 * qJD(3);
t96 = t20 * t71;
t22 = t45 * t32 - t43 * t74;
t94 = t22 * t71;
t93 = t23 * t71;
t92 = t25 * t71;
t30 = -pkin(7) + t54;
t91 = t30 * t47;
t31 = -t44 * pkin(4) - t42 * qJ(5) - pkin(3);
t90 = t31 * t71;
t89 = t71 * t42;
t88 = t71 * t44;
t85 = t47 * t42;
t84 = t47 * t44;
t78 = qJD(4) * t42;
t83 = t22 * t78 - t23 * t88;
t82 = t23 - t25;
t40 = t42 ^ 2;
t41 = t44 ^ 2;
t81 = t40 - t41;
t80 = t40 + t41;
t70 = t42 * t105 * t44;
t69 = 2 * t72;
t5 = t14 - t92;
t68 = -t5 - t98;
t16 = -t22 + t97;
t67 = t16 + t97;
t66 = -t14 - t98;
t6 = -t22 - t90;
t65 = t6 - t90;
t60 = t77 * t89;
t19 = -t31 + t99;
t59 = -t19 * t71 - t20 - t6;
t56 = t42 * t8 + t44 * t9;
t7 = t21 - t25;
t53 = t7 * t71 + t5 - t91;
t52 = t91 - t100;
t51 = qJD(4) * (-(pkin(3) + t99) * t71 - t16 - t20);
t50 = (-t47 - t105) * t43;
t48 = qJD(1) ^ 2;
t26 = t55 * t71;
t24 = t81 * t103;
t4 = t42 * t101 + t44 * t50;
t3 = -t44 * t101 + t42 * t50;
t12 = [0, 0, 0, 0, t69, qJ(2) * t69, 0, t100, t13 + t96, 0.2e1 * t60, -0.2e1 * t24, -t84, t85, 0, t42 * t51 - t52 * t44, t52 * t42 + t44 * t51, t53 * t44 + t59 * t78, -t80 * t96 - t49, t53 * t42 - t59 * t77, t5 * t19 + t56 * t20 + t49 * t30 + t6 * t7; 0, 0, 0, 0, -t48, -t48 * qJ(2), 0, -t43 * t105, -t102, 0, 0, 0, 0, 0, t4, -t3, t4, t80 * t102, t3, (-t71 * t56 - t5) * t45 + (-t71 * t6 + t49) * t43; 0, 0, 0, 0, 0, 0, 0, -t14 - t93, -t13 - t94, -0.2e1 * t60, 0.2e1 * t24, t84, -t85, 0, t66 * t44 + t67 * t78 + t83, (-t66 + t93) * t42 + (t22 + t67) * t77, t65 * t78 + (t68 + t92) * t44 + t83, t80 * t94 + t49, (-t22 - t65) * t77 + (-t71 * t82 + t68) * t42, t49 * pkin(7) - t56 * t22 + t5 * t31 - t82 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t81 * t105, 0, 0, 0, t16 * t89 - t10, t16 * t88 - t11, -t10 - (-t26 * t44 - t42 * t6) * t71, -((t9 - t73) * t42 + (t62 - t8) * t44) * t71, 0.2e1 * qJD(4) * qJD(5) + t11 - (-t26 * t42 + t44 * t6) * t71, -t2 * pkin(4) + t1 * qJ(5) + t9 * qJD(5) + t57 * t17 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, -t105 * t40 - t47, -t9 * qJD(4) - t6 * t89 + t2;];
tauc_reg = t12;
