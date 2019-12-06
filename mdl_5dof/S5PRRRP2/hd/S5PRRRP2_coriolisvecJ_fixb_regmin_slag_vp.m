% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:57
% EndTime: 2019-12-05 16:42:00
% DurationCPUTime: 0.39s
% Computational Cost: add. (510->104), mult. (939->139), div. (0->0), fcn. (434->4), ass. (0->77)
t43 = cos(qJ(3));
t76 = pkin(2) * qJD(3);
t59 = qJD(2) * t76;
t93 = qJD(1) * qJD(4) + t43 * t59;
t42 = cos(qJ(4));
t37 = qJD(2) + qJD(3);
t41 = sin(qJ(3));
t77 = pkin(2) * qJD(2);
t65 = t41 * t77;
t23 = t37 * pkin(7) + t65;
t40 = sin(qJ(4));
t84 = t40 * t23;
t12 = t42 * qJD(1) - t84;
t92 = qJD(5) - t12;
t80 = t93 * t42;
t2 = (qJD(5) - t84) * qJD(4) + t80;
t74 = qJD(4) * t42;
t5 = t23 * t74 + t93 * t40;
t91 = t2 * t42 + t5 * t40;
t38 = t40 ^ 2;
t39 = t42 ^ 2;
t90 = t37 * (t38 + t39);
t9 = -qJD(4) * pkin(4) + t92;
t13 = t40 * qJD(1) + t42 * t23;
t69 = qJD(4) * qJ(5);
t10 = t13 + t69;
t44 = qJD(4) ^ 2;
t89 = pkin(7) * t44;
t88 = t43 * pkin(2);
t25 = -t42 * pkin(4) - t40 * qJ(5) - pkin(3);
t87 = t25 * t37;
t86 = t37 * t40;
t85 = t37 * t42;
t83 = t44 * t40;
t35 = t44 * t42;
t64 = t43 * t77;
t24 = -t37 * pkin(3) - t64;
t56 = t41 * t59;
t82 = t24 * t74 + t40 * t56;
t73 = qJD(4) * t43;
t61 = t40 * t73;
t81 = t61 * t77 + t65 * t85;
t79 = t38 - t39;
t75 = qJD(4) * t40;
t72 = t10 * qJD(4);
t71 = t40 * qJD(5);
t70 = -qJD(2) - t37;
t36 = t37 ^ 2;
t67 = t40 * t36 * t42;
t66 = t43 * t76;
t63 = t37 * t75;
t62 = t37 * t74;
t51 = pkin(4) * t40 - qJ(5) * t42;
t4 = t56 + (t51 * qJD(4) - t71) * t37;
t60 = -t4 - t89;
t58 = t12 + t84;
t54 = (-qJD(3) + t37) * t77;
t53 = t70 * t76;
t52 = t10 * t42 + t40 * t9;
t50 = t13 * qJD(4) - t5;
t15 = pkin(4) * t75 - t42 * t69 - t71;
t11 = t41 * t76 + t15;
t32 = t41 * pkin(2) + pkin(7);
t49 = -t11 * t37 - t32 * t44 - t4;
t20 = t25 - t88;
t48 = t20 * t37 - t66;
t47 = -t40 * t72 + t9 * t74 + t91;
t46 = -t41 * t86 + t42 * t73;
t45 = (-t10 * t40 + t42 * t9) * qJD(4) + t91;
t33 = -pkin(3) - t88;
t22 = 0.2e1 * t40 * t62;
t17 = t24 * t75;
t16 = t51 * t37;
t14 = -0.2e1 * t79 * t37 * qJD(4);
t8 = -t64 + t87;
t6 = t8 * t75;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t35, -t83, 0, t35, t52 * qJD(4) + t2 * t40 - t5 * t42; 0, 0, 0, 0, 0, t41 * t53, t43 * t53, t22, t14, t35, -t83, 0, t33 * t63 - t32 * t35 + t17 + (t70 * t42 * t41 - t61) * t76, t32 * t83 + t33 * t62 - t46 * t76 + t82, t49 * t42 + t48 * t75 + t6, t66 * t90 + t47, t49 * t40 + (-t48 - t8) * t74, t8 * t11 + t4 * t20 + t45 * t32 + t52 * t66; 0, 0, 0, 0, 0, t41 * t54, t43 * t54, t22, t14, t35, -t83, 0, -pkin(3) * t63 + t17 + (-t56 - t89) * t42 + t81, -pkin(3) * t62 + pkin(7) * t83 + t46 * t77 + t82, t25 * t63 + t6 + (-t15 * t37 + t60) * t42 + t81, -t64 * t90 + t47, (-t64 - t8 - t87) * t74 + ((-t15 + t65) * t37 + t60) * t40, t8 * t15 + t4 * t25 + (-t41 * t8 - t52 * t43) * t77 + t45 * pkin(7); 0, 0, 0, 0, 0, 0, 0, -t67, t79 * t36, 0, 0, 0, -t24 * t86 + t50, t58 * qJD(4) - t24 * t85 - t80, (t16 * t42 - t40 * t8) * t37 + t50, 0, (t16 * t40 + t42 * t8) * t37 + (0.2e1 * qJD(5) - t58) * qJD(4) + t80, -t5 * pkin(4) + t2 * qJ(5) + t92 * t10 - t9 * t13 - t8 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, -t38 * t36 - t44, t8 * t86 + t5 - t72;];
tauc_reg = t1;
