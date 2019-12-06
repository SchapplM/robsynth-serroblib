% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:18
% EndTime: 2019-12-05 15:11:20
% DurationCPUTime: 0.49s
% Computational Cost: add. (306->79), mult. (956->122), div. (0->0), fcn. (872->6), ass. (0->72)
t53 = sin(qJ(4));
t76 = t53 * qJD(5);
t49 = t53 ^ 2;
t55 = cos(qJ(4));
t50 = t55 ^ 2;
t88 = t49 + t50;
t96 = t88 * qJD(3);
t97 = pkin(6) * t96 + t76;
t87 = t55 * qJ(5);
t94 = t53 * pkin(4);
t40 = -t87 + t94;
t95 = t40 / 0.2e1;
t52 = cos(pkin(8));
t51 = sin(pkin(8));
t56 = cos(qJ(3));
t90 = t51 * t56;
t27 = t52 * t55 + t53 * t90;
t93 = t27 * t53;
t28 = -t52 * t53 + t55 * t90;
t92 = t28 * t55;
t54 = sin(qJ(3));
t91 = t51 * t54;
t89 = t56 * t54;
t86 = qJD(3) * t53;
t85 = qJD(3) * t54;
t84 = qJD(3) * t55;
t83 = qJD(3) * t56;
t82 = qJD(4) * t53;
t48 = qJD(4) * t55;
t81 = qJD(4) * t56;
t64 = -t55 * pkin(4) - t53 * qJ(5);
t37 = -pkin(3) + t64;
t13 = t37 * t55 + t40 * t53;
t80 = t13 * qJD(3);
t14 = -t37 * t53 + t40 * t55;
t79 = t14 * qJD(3);
t46 = t50 - t49;
t78 = t46 * qJD(3);
t77 = t49 * qJD(3);
t75 = t55 * qJD(5);
t74 = qJD(4) * qJ(5);
t73 = pkin(3) * t86;
t72 = pkin(3) * t84;
t71 = pkin(6) * t82;
t70 = pkin(6) * t48;
t69 = t53 * t85;
t68 = t54 * t84;
t67 = t37 * t86;
t4 = (-t92 / 0.2e1 - t93 / 0.2e1 + t90 / 0.2e1) * t56 + (t50 / 0.2e1 - 0.1e1 / 0.2e1 + t49 / 0.2e1) * t51 * t54 ^ 2;
t6 = -t51 ^ 2 * t89 + (t92 + t93) * t91;
t63 = -t6 * qJD(1) - t4 * qJD(2);
t15 = (-0.1e1 + t88) * t89;
t62 = -t4 * qJD(1) + t15 * qJD(2);
t61 = t87 / 0.2e1 - t94 / 0.2e1;
t60 = t27 * qJD(4) + t51 * t68;
t12 = -t28 * qJD(4) + t51 * t69;
t29 = t54 * t82 - t55 * t83;
t30 = t54 * t48 + t53 * t83;
t59 = t95 + t61;
t10 = t59 * t56;
t2 = t59 * t91;
t58 = t37 * t40 * qJD(3) + t2 * qJD(1) - t10 * qJD(2);
t57 = t64 * qJD(4) + t75;
t47 = t53 * t84;
t32 = -t53 * t81 - t68;
t31 = -t55 * t81 + t69;
t23 = t29 * t51;
t22 = t30 * t51;
t11 = (-t40 / 0.2e1 + t61) * t56;
t3 = (-t61 + t95) * t91;
t1 = t4 * qJD(3);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, -t51 * t83, t51 * t85, 0, 0, 0, 0, 0, t23, t22, t23, -t96 * t91, -t22, t3 * qJD(4) + (t37 * t83 - t97 * t54) * t51 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t60, t12, 0, -t60, t3 * qJD(3) + (-t28 * pkin(4) - t27 * qJ(5)) * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * qJD(3); 0, 0, 0, -t85, -t83, 0, 0, 0, 0, 0, t32, t31, t32, t56 * t96, -t31, t11 * qJD(4) + t37 * t85 + t97 * t56 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, -t30, 0, -t29, t11 * qJD(3) + t54 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(4) - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * qJD(4) - t62; 0, 0, 0, 0, 0, t53 * t48, t46 * qJD(4), 0, 0, 0, -pkin(3) * t82, -pkin(3) * t48, -t14 * qJD(4) + t53 * t75, 0, -t13 * qJD(4) + t49 * qJD(5), (qJD(4) * t40 - t76) * t37; 0, 0, 0, 0, 0, t47, t78, t48, -t82, 0, -t70 - t73, t71 - t72, -t70 - t79, t57, -t71 - t80, pkin(6) * t57 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t48, t77, -t67 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(3); 0, 0, 0, 0, 0, -t47, -t78, 0, 0, 0, t73, t72, t79, 0, t80, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, -t77, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
