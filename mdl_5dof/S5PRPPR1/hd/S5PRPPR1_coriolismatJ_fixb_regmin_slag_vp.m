% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:12
% EndTime: 2019-12-05 15:22:15
% DurationCPUTime: 0.48s
% Computational Cost: add. (338->78), mult. (990->132), div. (0->0), fcn. (981->6), ass. (0->84)
t59 = sin(pkin(9));
t63 = sin(qJ(5));
t100 = t63 * t59;
t61 = cos(pkin(9));
t64 = cos(qJ(5));
t97 = t64 * t61;
t66 = -t97 + t100;
t104 = t61 ^ 2;
t105 = t59 ^ 2;
t106 = t104 + t105;
t62 = cos(pkin(8));
t58 = t62 ^ 2;
t60 = sin(pkin(8));
t102 = pkin(6) * t60;
t101 = t60 * t62;
t99 = t63 * t61;
t98 = t64 * t59;
t32 = t66 * t60;
t82 = t32 * qJD(2);
t25 = t62 * t82;
t29 = t32 * qJD(5);
t96 = -t25 + t29;
t42 = -pkin(3) * t62 - qJ(4) * t60 - pkin(2);
t95 = qJ(3) * t62;
t27 = t59 * t42 + t61 * t95;
t57 = t60 ^ 2;
t51 = t57 + t58;
t38 = t61 * t42;
t17 = -t61 * t102 + t38 + (-qJ(3) * t59 - pkin(4)) * t62;
t19 = -t102 * t59 + t27;
t3 = -t17 * t64 + t19 * t63;
t67 = -t98 - t99;
t30 = t67 * t60;
t68 = (pkin(4) * t59 + qJ(3)) * t60;
t1 = -t3 * t62 + t30 * t68;
t94 = t1 * qJD(2);
t4 = -t17 * t63 - t19 * t64;
t2 = t32 * t68 + t4 * t62;
t93 = t2 * qJD(2);
t26 = -t59 * t95 + t38;
t6 = (t26 * t61 + t27 * t59) * t60;
t92 = t6 * qJD(2);
t7 = t30 ^ 2 - t32 ^ 2;
t91 = t7 * qJD(2);
t10 = t60 * t30 + t58 * t67;
t90 = t10 * qJD(2);
t11 = -t60 * t32 - t58 * t66;
t89 = t11 * qJD(2);
t13 = t67 * t62;
t88 = t13 * qJD(2);
t14 = t66 * t62;
t87 = t14 * qJD(2);
t69 = t105 / 0.2e1 + t104 / 0.2e1;
t21 = (-0.1e1 / 0.2e1 + t69) * t101;
t86 = t21 * qJD(2);
t85 = t30 * qJD(2);
t84 = t30 * qJD(4);
t83 = t30 * qJD(5);
t81 = t32 * qJD(4);
t33 = (0.1e1 / 0.2e1 + t69) * t60;
t80 = t33 * qJD(2);
t39 = t106 * t57;
t79 = t39 * qJD(2);
t40 = t51 * t59;
t78 = t40 * qJD(2);
t41 = t51 * t61;
t77 = t41 * qJD(2);
t49 = t51 * qJ(3);
t76 = t49 * qJD(2);
t75 = t51 * qJD(2);
t74 = qJD(2) * t101;
t73 = qJD(4) * t101;
t72 = t30 * t82;
t71 = t59 * t74;
t70 = t61 * t74;
t5 = t57 * qJ(3) + (-t26 * t59 + t27 * t61) * t62;
t65 = -qJD(1) * t21 - qJD(2) * t5;
t34 = (0.1e1 / 0.2e1 - t106 / 0.2e1) * t60;
t24 = t62 * t85;
t20 = t21 * qJD(3);
t16 = -t13 / 0.2e1 + (-t99 / 0.2e1 - t98 / 0.2e1) * t62;
t15 = (-t66 / 0.2e1 - t97 / 0.2e1 + t100 / 0.2e1) * t62;
t12 = -t24 + t83;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t51 * qJD(3), t49 * qJD(3), qJD(3) * t40 + t61 * t73, qJD(3) * t41 - t59 * t73, t39 * qJD(4), qJD(3) * t5 - qJD(4) * t6, -t30 * t29, t7 * qJD(5), -t62 * t83, -t62 * t29, 0, -qJD(3) * t10 - qJD(5) * t2 - t62 * t81, qJD(3) * t11 + qJD(5) * t1 + t62 * t84; 0, 0, 0, 0, 0, 0, t75, t76, t78, t77, 0, qJD(4) * t34 - t65, 0, 0, 0, 0, 0, qJD(5) * t16 - t90, qJD(5) * t15 + t89; 0, 0, 0, 0, 0, 0, 0, 0, t70, -t71, t79, qJD(3) * t34 - t92, 0, 0, 0, 0, 0, -t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t91, t12, t96, 0, qJD(3) * t16 + qJD(5) * t4 - t93, qJD(3) * t15 + qJD(5) * t3 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t75, -t76, -t78, -t77, 0, -qJD(4) * t33 + t65, 0, 0, 0, 0, 0, -qJD(5) * t13 + t90, -qJD(5) * t14 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t67 - t88, qJD(5) * t66 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t70, t71, -t79, qJD(3) * t33 + t92, 0, 0, 0, 0, 0, -t96, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t91, t24, t25, 0, qJD(3) * t13 + t81 + t93, qJD(3) * t14 - t84 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
