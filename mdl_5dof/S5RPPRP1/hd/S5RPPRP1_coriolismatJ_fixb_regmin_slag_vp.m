% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:52
% EndTime: 2020-01-03 11:25:55
% DurationCPUTime: 0.51s
% Computational Cost: add. (761->75), mult. (1453->128), div. (0->0), fcn. (1262->6), ass. (0->83)
t53 = cos(pkin(8));
t52 = sin(pkin(8));
t31 = -cos(pkin(7)) * pkin(1) - pkin(2) - t52 * pkin(6) - t53 * pkin(3);
t55 = cos(qJ(4));
t27 = t55 * t31;
t95 = qJ(5) * t52;
t60 = -t55 * t95 + t27;
t46 = sin(pkin(7)) * pkin(1) + qJ(3);
t54 = sin(qJ(4));
t99 = t46 * t54;
t20 = (-pkin(4) - t99) * t53 + t60;
t74 = t53 * t99;
t21 = t60 - t74;
t104 = -t20 + t21;
t50 = t54 ^ 2;
t51 = t55 ^ 2;
t103 = t50 + t51;
t102 = t21 / 0.2e1;
t100 = t20 * t54;
t48 = t52 ^ 2;
t98 = t48 * t55;
t96 = t52 * t55;
t41 = t53 ^ 2 + t48;
t66 = -t20 / 0.2e1 + t102;
t73 = -t53 * pkin(4) / 0.2e1;
t2 = (t73 + t66) * t96;
t94 = t2 * qJD(1);
t4 = t104 * t52 * t54;
t93 = t4 * qJD(1);
t5 = (t73 - t66) * t54;
t92 = t5 * qJD(1);
t24 = -t55 * t53 * t46 - t54 * t31;
t22 = -t54 * t95 - t24;
t8 = (t20 * t55 + t22 * t54) * t52;
t91 = t8 * qJD(1);
t90 = qJD(1) * t53;
t89 = qJD(3) * t53;
t88 = qJD(4) * t54;
t87 = qJD(4) * t55;
t23 = -t27 + t74;
t11 = -t23 * t53 - t48 * t99;
t86 = t11 * qJD(1);
t12 = t24 * t53 - t46 * t98;
t85 = t12 * qJD(1);
t65 = t50 / 0.2e1 + t51 / 0.2e1;
t25 = (-0.1e1 / 0.2e1 + t65) * t53 * t52;
t84 = t25 * qJD(1);
t26 = t41 * t46;
t83 = t26 * qJD(1);
t28 = (0.1e1 / 0.2e1 + t65) * t52;
t82 = t28 * qJD(1);
t32 = t103 * t48;
t81 = t32 * qJD(1);
t33 = (t50 - t51) * t48;
t80 = t33 * qJD(1);
t34 = t41 * t54;
t79 = t34 * qJD(1);
t35 = t41 * t55;
t78 = t35 * qJD(1);
t77 = t41 * qJD(1);
t76 = pkin(4) * t96;
t75 = t54 * t98;
t72 = t52 * t88;
t71 = t53 * t88;
t70 = t52 * t87;
t69 = t53 * t87;
t68 = t54 * t90;
t67 = t55 * t90;
t64 = qJD(1) * t76;
t63 = pkin(4) * t70;
t62 = qJD(1) * t75;
t61 = -qJD(4) + t90;
t30 = (pkin(4) * t54 + t46) * t52;
t3 = t104 * t22 + t30 * t76;
t59 = t3 * qJD(1) + t2 * qJD(2);
t7 = t30 * t52 + (t22 * t55 - t100) * t53;
t58 = -t7 * qJD(1) - t25 * qJD(2);
t57 = t61 * t54;
t56 = t61 * t55;
t29 = (0.1e1 / 0.2e1 - t103 / 0.2e1) * t52;
t6 = -t100 / 0.2e1 + (t102 + t73) * t54;
t1 = t25 * qJD(3) + t2 * qJD(4);
t9 = [0, 0, 0, 0, 0, 0, t41 * qJD(3), t26 * qJD(3), -qJD(4) * t75, t33 * qJD(4), t52 * t71, t52 * t69, 0, t34 * qJD(3) - t12 * qJD(4), t35 * qJD(3) + t11 * qJD(4), -t4 * qJD(4) + t32 * qJD(5), t7 * qJD(3) + t3 * qJD(4) - t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, t77, t83, 0, 0, 0, 0, 0, t79, t78, 0, t6 * qJD(4) + t29 * qJD(5) - t58; 0, 0, 0, 0, 0, 0, 0, 0, -t62, t80, t52 * t57, t52 * t56, 0, t24 * qJD(4) - t85, t23 * qJD(4) + t86, pkin(4) * t72 - t93, -t22 * pkin(4) * qJD(4) + t6 * qJD(3) + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t29 * qJD(3) - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t72, 0, -t63 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t77, -t83, 0, 0, 0, 0, 0, t71 - t79, t69 - t78, 0, -t5 * qJD(4) - t28 * qJD(5) + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56, 0, -pkin(4) * t88 - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, t62, -t80, -t52 * t68, -t52 * t67, 0, -t54 * t89 + t85, -t55 * t89 - t86, t93, t5 * qJD(3) - qJD(5) * t76 - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, 0, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t28 * qJD(3) + t63 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
