% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:45
% EndTime: 2022-01-20 09:51:47
% DurationCPUTime: 0.58s
% Computational Cost: add. (611->102), mult. (1359->135), div. (0->0), fcn. (1245->8), ass. (0->92)
t67 = cos(pkin(8));
t69 = sin(qJ(2));
t106 = t67 * t69;
t70 = cos(qJ(2));
t110 = t70 * pkin(1);
t60 = pkin(2) + t110;
t65 = sin(pkin(8));
t78 = pkin(1) * t106 + t65 * t60;
t36 = qJ(4) + t78;
t58 = t65 * pkin(2) + qJ(4);
t64 = sin(pkin(9));
t62 = t64 ^ 2;
t66 = cos(pkin(9));
t63 = t66 ^ 2;
t123 = (t63 / 0.2e1 + t62 / 0.2e1) * (t36 + t58);
t68 = sin(qJ(5));
t104 = t68 * t66;
t109 = cos(qJ(5));
t83 = t109 * t64;
t48 = t83 + t104;
t111 = t66 * pkin(4);
t57 = t65 * t69 * pkin(1);
t81 = t67 * t60 - t57;
t77 = -pkin(3) - t81;
t28 = t77 - t111;
t84 = -t67 * pkin(2) - pkin(3);
t51 = t84 - t111;
t86 = t51 / 0.2e1 + t28 / 0.2e1;
t121 = t86 * t48;
t105 = t68 * t64;
t82 = t109 * t66;
t46 = -t82 + t105;
t14 = t46 ^ 2 - t48 ^ 2;
t91 = qJD(1) + qJD(2);
t120 = t91 * t14;
t119 = t91 * t46;
t79 = t91 * t48;
t56 = t62 + t63;
t118 = t91 * t56;
t107 = t65 * t70;
t44 = (t106 + t107) * pkin(1);
t108 = t44 * t66;
t45 = t67 * t110 - t57;
t16 = t56 * t45;
t52 = t56 * qJD(4);
t103 = t16 * qJD(2) + t52;
t102 = pkin(1) * qJD(1);
t101 = pkin(1) * qJD(2);
t1 = t36 * t16 + t77 * t44;
t100 = t1 * qJD(1);
t10 = -t81 * t44 + t78 * t45;
t99 = t10 * qJD(1);
t12 = t56 * t36;
t98 = t12 * qJD(1);
t97 = t16 * qJD(1);
t96 = t46 * qJD(1);
t95 = t46 * qJD(2);
t40 = t46 * qJD(5);
t94 = t48 * qJD(1);
t93 = t48 * qJD(2);
t92 = t48 * qJD(5);
t90 = t28 * t96;
t89 = t28 * t94;
t88 = t44 * t96;
t87 = t44 * t94;
t33 = t56 * t58;
t80 = pkin(1) * t91;
t73 = (t107 / 0.2e1 + t106 / 0.2e1) * pkin(1);
t8 = t73 - t123;
t76 = t8 * qJD(1) - t33 * qJD(2);
t72 = (-t104 / 0.2e1 - t83 / 0.2e1) * t45;
t2 = t72 - t121;
t75 = -t2 * qJD(1) + t51 * t93;
t71 = (-t82 / 0.2e1 + t105 / 0.2e1) * t45;
t3 = t86 * t46 + t71;
t74 = -t3 * qJD(1) - t51 * t95;
t61 = t66 * pkin(7);
t43 = t48 * qJD(4);
t39 = t46 * qJD(4);
t38 = t66 * t58 + t61;
t37 = (-pkin(7) - t58) * t64;
t23 = t66 * t36 + t61;
t22 = (-pkin(7) - t36) * t64;
t21 = t46 * t92;
t20 = t44 * t93;
t19 = t44 * t95;
t13 = t14 * qJD(5);
t11 = t46 * t79;
t9 = t73 + t123;
t5 = t72 + t121;
t4 = t71 - (t28 + t51) * t46 / 0.2e1;
t6 = [0, 0, 0, 0, -t69 * t101, -t70 * t101, t10 * qJD(2), -qJD(2) * t108, t103, qJD(2) * t1 + qJD(4) * t12, -t21, t13, 0, 0, 0, t28 * t92 + t19, -t28 * t40 + t20; 0, 0, 0, 0, -t69 * t80, -t70 * t80, t99 + (-t44 * t67 + t45 * t65) * qJD(2) * pkin(2), -t91 * t108, t97 + t103, t100 + (t33 * t45 + t44 * t84) * qJD(2) + t9 * qJD(4), -t21, t13, 0, 0, 0, qJD(5) * t5 + t19 + t88, qJD(5) * t4 + t20 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t118, qJD(2) * t9 + t98, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t120, -t40, -t92, 0, t89 + t5 * qJD(2) + (-t109 * t23 - t22 * t68) * qJD(5), -t90 + t4 * qJD(2) + (-t109 * t22 + t23 * t68) * qJD(5); 0, 0, 0, 0, t69 * t102, t70 * t102, -t99, qJD(1) * t108, t52 - t97, -qJD(4) * t8 - t100, -t21, t13, 0, 0, 0, -qJD(5) * t2 - t88, -qJD(5) * t3 - t87; 0, 0, 0, 0, 0, 0, 0, 0, t52, t33 * qJD(4), -t21, t13, 0, 0, 0, t51 * t92, -t51 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t118, -t76, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t120, -t40, -t92, 0, (-t109 * t38 - t37 * t68) * qJD(5) + t75, (-t109 * t37 + t38 * t68) * qJD(5) + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t40; 0, 0, 0, 0, 0, 0, 0, 0, -t118, qJD(2) * t8 - t98, 0, 0, 0, 0, 0, t92, -t40; 0, 0, 0, 0, 0, 0, 0, 0, -t118, t76, 0, 0, 0, 0, 0, t92, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t120, 0, 0, 0, qJD(2) * t2 - t43 - t89, qJD(2) * t3 + t39 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t120, 0, 0, 0, -t43 - t75, t39 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
