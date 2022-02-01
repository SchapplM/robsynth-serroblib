% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:47
% EndTime: 2022-01-20 09:12:49
% DurationCPUTime: 0.49s
% Computational Cost: add. (406->77), mult. (1040->131), div. (0->0), fcn. (1040->8), ass. (0->80)
t59 = sin(pkin(9));
t61 = cos(pkin(9));
t69 = t59 ^ 2 / 0.2e1 + t61 ^ 2 / 0.2e1;
t105 = -0.1e1 / 0.2e1 + t69;
t64 = cos(qJ(5));
t94 = t64 * t61;
t63 = sin(qJ(5));
t97 = t63 * t59;
t104 = t94 - t97;
t62 = cos(pkin(8));
t58 = t62 ^ 2;
t60 = sin(pkin(8));
t100 = pkin(6) * t60;
t55 = sin(pkin(7)) * pkin(1) + qJ(3);
t99 = t55 * t62;
t98 = t62 * t60;
t96 = t63 * t61;
t95 = t64 * t59;
t33 = t104 * t60;
t79 = t33 * qJD(1);
t27 = t62 * t79;
t29 = t33 * qJD(5);
t93 = t27 - t29;
t40 = -cos(pkin(7)) * pkin(1) - pkin(2) - t60 * qJ(4) - t62 * pkin(3);
t21 = t59 * t40 + t61 * t99;
t57 = t60 ^ 2;
t51 = t57 + t58;
t35 = t61 * t40;
t13 = -t61 * t100 + t35 + (-t55 * t59 - pkin(4)) * t62;
t14 = -t100 * t59 + t21;
t3 = -t13 * t64 + t14 * t63;
t67 = -t95 - t96;
t31 = t67 * t60;
t68 = (pkin(4) * t59 + t55) * t60;
t1 = -t3 * t62 + t31 * t68;
t92 = qJD(1) * t1;
t4 = -t13 * t63 - t14 * t64;
t2 = -t33 * t68 + t4 * t62;
t91 = qJD(1) * t2;
t20 = -t59 * t99 + t35;
t6 = (t20 * t61 + t21 * t59) * t60;
t90 = qJD(1) * t6;
t7 = t31 ^ 2 - t33 ^ 2;
t89 = t7 * qJD(1);
t10 = t60 * t31 + t58 * t67;
t88 = qJD(1) * t10;
t11 = t104 * t58 + t60 * t33;
t87 = qJD(1) * t11;
t86 = qJD(4) * t62;
t15 = t67 * t62;
t85 = t15 * qJD(1);
t16 = t104 * t62;
t84 = t16 * qJD(1);
t23 = t105 * t98;
t83 = t23 * qJD(1);
t30 = t51 * t55;
t82 = t30 * qJD(1);
t81 = t31 * qJD(1);
t80 = t31 * qJD(5);
t36 = (0.1e1 / 0.2e1 + t69) * t60;
t78 = t36 * qJD(1);
t41 = t51 * t59;
t77 = t41 * qJD(1);
t42 = t51 * t61;
t76 = t42 * qJD(1);
t75 = t51 * qJD(1);
t74 = qJD(1) * t98;
t73 = t60 * t86;
t72 = t31 * t79;
t71 = t59 * t74;
t70 = t61 * t74;
t5 = t55 * t57 + (-t20 * t59 + t21 * t61) * t62;
t65 = -qJD(1) * t5 - qJD(2) * t23;
t37 = t105 * t60;
t26 = t62 * t81;
t22 = t23 * qJD(3);
t18 = -t15 / 0.2e1 + (-t96 / 0.2e1 - t95 / 0.2e1) * t62;
t17 = (t104 / 0.2e1 - t94 / 0.2e1 + t97 / 0.2e1) * t62;
t12 = -t26 + t80;
t8 = [0, 0, 0, 0, 0, t51 * qJD(3), t30 * qJD(3), qJD(3) * t41 + t61 * t73, qJD(3) * t42 - t59 * t73, qJD(3) * t5 - qJD(4) * t6, t31 * t29, t7 * qJD(5), -t62 * t80, t62 * t29, 0, -qJD(3) * t10 - qJD(5) * t2 + t33 * t86, qJD(3) * t11 + qJD(5) * t1 + t31 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t75, t82, t77, t76, -qJD(4) * t37 - t65, 0, 0, 0, 0, 0, qJD(5) * t18 - t88, qJD(5) * t17 + t87; 0, 0, 0, 0, 0, 0, 0, t70, -t71, -qJD(3) * t37 - t90, 0, 0, 0, 0, 0, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t89, t12, t93, 0, qJD(3) * t18 + qJD(5) * t4 - t91, qJD(3) * t17 + qJD(5) * t3 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t80; 0, 0, 0, 0, 0, -t75, -t82, -t77, -t76, -qJD(4) * t36 + t65, 0, 0, 0, 0, 0, -qJD(5) * t15 + t88, qJD(5) * t16 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t67 - t85, -qJD(5) * t104 + t84; 0, 0, 0, 0, 0, 0, 0, -t70, t71, qJD(3) * t36 + t90, 0, 0, 0, 0, 0, -t93, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t89, t26, -t27, 0, qJD(3) * t15 - qJD(4) * t33 + t91, -qJD(3) * t16 - qJD(4) * t31 - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
