% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:18
% DurationCPUTime: 0.41s
% Computational Cost: add. (408->92), mult. (1059->136), div. (0->0), fcn. (705->6), ass. (0->75)
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t38 = sin(pkin(8));
t81 = qJD(1) * t38;
t27 = t41 * qJD(2) + t43 * t81;
t98 = t27 * qJD(3);
t21 = qJD(3) * pkin(6) + t27;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t39 = cos(pkin(8));
t80 = qJD(1) * t39;
t50 = t40 * t21 + t42 * t80;
t97 = qJD(5) + t50;
t5 = -qJD(4) * pkin(4) + t97;
t93 = t43 * qJD(2) - t41 * t81;
t23 = t93 * qJD(3);
t69 = qJD(3) * qJD(4);
t95 = -0.2e1 * t69;
t44 = qJD(4) ^ 2;
t45 = qJD(3) ^ 2;
t94 = (t44 + t45) * t41;
t63 = t40 * t80;
t10 = t42 * t21 - t63;
t6 = qJD(4) * qJ(5) + t10;
t92 = pkin(6) * t44;
t91 = t38 * t41;
t90 = t38 * t43;
t89 = t45 * t41;
t88 = t45 * t43;
t75 = qJD(4) * t40;
t87 = t42 * t98 + t75 * t93;
t52 = pkin(4) * t40 - qJ(5) * t42;
t19 = t52 * qJD(4) - t40 * qJD(5);
t86 = t19 - t27;
t36 = t40 ^ 2;
t37 = t42 ^ 2;
t85 = t36 - t37;
t84 = t36 + t37;
t82 = qJD(3) * pkin(3);
t30 = -t42 * pkin(4) - t40 * qJ(5) - pkin(3);
t77 = qJD(3) * t30;
t11 = -t93 + t77;
t79 = qJD(3) * t11;
t78 = qJD(3) * t19;
t76 = qJD(3) * t40;
t74 = qJD(4) * t42;
t68 = t42 * t90;
t67 = t40 * t45 * t42;
t66 = t38 * t88;
t65 = qJD(4) * t63 - t21 * t74 - t40 * t23;
t64 = qJD(3) * t91;
t4 = t78 + t98;
t61 = -t4 - t92;
t60 = -t98 - t92;
t59 = t40 * t69;
t20 = -t93 - t82;
t58 = t20 - t82;
t57 = t11 + t77;
t55 = t42 * t64;
t54 = t43 * t95;
t53 = t40 * t5 + t42 * t6;
t51 = t10 * qJD(4) + t65;
t24 = t39 * t42 + t40 * t90;
t7 = -qJD(4) * t68 + (qJD(4) * t39 + t64) * t40;
t47 = t7 * qJD(4) - t42 * t66 + t59 * t91;
t16 = t42 * t23;
t2 = t16 + (qJD(5) - t50) * qJD(4);
t46 = t2 * t42 - t65 * t40 + (-t40 * t6 + t42 * t5) * qJD(4);
t28 = t52 * qJD(3);
t25 = -t39 * t40 + t68;
t13 = t40 * t54 - t42 * t94;
t12 = t40 * t94 + t42 * t54;
t8 = -qJD(4) * t24 - t55;
t1 = -t40 * t66 + (t8 - t55) * qJD(4);
t3 = [0, 0, 0, -t66, t38 * t89, 0, 0, 0, 0, 0, t47, -t1, t47, (-t40 * t7 + t42 * t8 + (t24 * t42 - t25 * t40) * qJD(4)) * qJD(3), t1, t2 * t25 - t65 * t24 - t5 * t7 + t6 * t8 + (t4 * t41 + t43 * t79) * t38; 0, 0, 0, -t89, -t88, 0, 0, 0, 0, 0, t13, t12, t13, t84 * t88, -t12, (t53 * qJD(3) - t4) * t43 + (t46 + t79) * t41; 0, 0, 0, 0, 0, 0.2e1 * t42 * t59, t85 * t95, t44 * t42, -t44 * t40, 0, t60 * t42 + t58 * t75 + t87, (-t60 - t98) * t40 + (t93 + t58) * t74, t57 * t75 + (t61 - t78) * t42 + t87, -t84 * t23 + t46, (-t93 - t57) * t74 + (-t86 * qJD(3) + t61) * t40, t46 * pkin(6) + t86 * t11 + t4 * t30 - t53 * t93; 0, 0, 0, 0, 0, -t67, t85 * t45, 0, 0, 0, -t20 * t76 + t51, -t20 * t42 * qJD(3) - t16, (-t11 * t40 + t28 * t42) * qJD(3) + t51, 0, t16 + (t11 * t42 + t28 * t40) * qJD(3) + 0.2e1 * qJD(4) * qJD(5), t65 * pkin(4) + t2 * qJ(5) - t5 * t10 - t11 * t28 + t97 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, -t36 * t45 - t44, -t6 * qJD(4) + t11 * t76 - t65;];
tauc_reg = t3;
