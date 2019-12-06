% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:53
% EndTime: 2019-12-05 17:04:55
% DurationCPUTime: 0.52s
% Computational Cost: add. (424->106), mult. (1119->153), div. (0->0), fcn. (834->6), ass. (0->93)
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t52 = -t58 ^ 2 + t61 ^ 2;
t84 = -qJD(3) - qJD(4);
t71 = qJD(2) - t84;
t115 = t71 * t52;
t62 = cos(qJ(4));
t63 = cos(qJ(3));
t105 = t62 * t63;
t59 = sin(qJ(4));
t60 = sin(qJ(3));
t108 = t59 * t60;
t40 = (t105 - t108) * pkin(2);
t111 = t62 * pkin(3);
t80 = -t111 / 0.2e1;
t114 = -t40 / 0.2e1 + t80;
t112 = t58 / 0.2e1;
t17 = t59 * pkin(3);
t110 = t63 * pkin(2);
t106 = t62 * t60;
t107 = t59 * t63;
t39 = (t106 + t107) * pkin(2);
t109 = t39 * t61;
t72 = pkin(3) + t110;
t66 = t59 * t72;
t36 = pkin(2) * t106 + t66;
t33 = t36 * qJD(4);
t37 = t39 * qJD(3);
t104 = -t37 - t33;
t103 = pkin(2) * qJD(2);
t102 = pkin(2) * qJD(3);
t68 = t17 / 0.2e1 + t36 / 0.2e1;
t64 = -t39 / 0.2e1 + t68;
t2 = t64 * t58;
t101 = qJD(2) * t2;
t3 = t64 * t61;
t100 = qJD(2) * t3;
t44 = t62 * t72;
t35 = pkin(2) * t108 - t44;
t22 = t35 * t112;
t67 = t80 + t40 / 0.2e1;
t5 = t58 * t67 + t22;
t99 = qJD(2) * t5;
t23 = t35 * t61 / 0.2e1;
t7 = t61 * t67 + t23;
t98 = qJD(2) * t7;
t97 = qJD(2) * t35;
t96 = qJD(2) * t36;
t95 = qJD(2) * t39;
t94 = qJD(2) * t40;
t93 = qJD(2) * t58;
t92 = qJD(2) * t61;
t91 = qJD(3) * t59;
t90 = qJD(4) * t59;
t89 = qJD(5) * t58;
t57 = qJD(5) * t61;
t88 = qJD(5) * t62;
t87 = t17 * qJD(2);
t18 = t44 / 0.2e1 + (-t110 / 0.2e1 + pkin(3) / 0.2e1) * t62;
t86 = t18 * qJD(2);
t85 = -qJD(2) - qJD(3);
t83 = pkin(3) * t91;
t82 = qJD(3) * t111;
t81 = pkin(3) * t90;
t79 = t35 * t93;
t78 = t35 * t92;
t77 = t36 * t93;
t76 = t39 * t93;
t75 = t58 * t91;
t74 = t58 * t88;
t73 = t61 * t88;
t70 = pkin(2) * t85;
t69 = t84 * t59;
t65 = t84 * t111;
t54 = pkin(6) + t17;
t53 = t58 * t57;
t51 = t58 * t81;
t43 = t52 * qJD(5);
t38 = t40 * qJD(3);
t34 = pkin(6) + t36;
t32 = t35 * qJD(4);
t31 = t58 * t37;
t26 = t58 * t33;
t19 = t71 * t61 * t58;
t14 = t80 - t44 / 0.2e1 + (t108 - t105 / 0.2e1) * pkin(2);
t13 = -t17 / 0.2e1 - t66 / 0.2e1 + (-t106 - t107 / 0.2e1) * pkin(2);
t12 = 0.2e1 * t23;
t11 = 0.2e1 * t22;
t8 = t114 * t61 + t23;
t6 = t114 * t58 + t22;
t4 = -t109 / 0.2e1 - t68 * t61;
t1 = t112 * t39 + t58 * t68;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t60 * t102, -t63 * t102, 0, t104, -t38 + t32, t53, t43, 0, 0, 0, t104 * t61 + t35 * t89, t35 * t57 + t26 + t31; 0, 0, 0, 0, 0, t60 * t70, t63 * t70, 0, qJD(4) * t13 - t37 - t95, qJD(4) * t14 - t38 - t94, t53, t43, 0, 0, 0, qJD(4) * t4 + qJD(5) * t6 + t109 * t85, qJD(4) * t1 + qJD(5) * t8 + t31 + t76; 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t13 - t33 - t96, qJD(3) * t14 + t32 + t97, t53, t43, 0, 0, 0, qJD(3) * t4 + qJD(5) * t11 + (-qJD(2) - qJD(4)) * t61 * t36, qJD(3) * t1 + qJD(5) * t12 + t26 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t115, t57, -t89, 0, qJD(3) * t6 + qJD(4) * t11 - t34 * t57 + t79, qJD(3) * t8 + qJD(4) * t12 + t34 * t89 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t60 * t103, t63 * t103, 0, -qJD(4) * t17 + t95, -qJD(4) * t18 + t94, t53, t43, 0, 0, 0, -qJD(4) * t3 + qJD(5) * t5 + t39 * t92, qJD(4) * t2 + qJD(5) * t7 - t76; 0, 0, 0, 0, 0, 0, 0, 0, -t81, -qJD(4) * t111, t53, t43, 0, 0, 0, (-t61 * t90 - t74) * pkin(3), -pkin(3) * t73 + t51; 0, 0, 0, 0, 0, 0, 0, 0, pkin(3) * t69 - t87, t65 - t86, t53, t43, 0, 0, 0, -t100 + (t61 * t69 - t74) * pkin(3), t101 + t51 + (-t73 + t75) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t115, t57, -t89, 0, -t54 * t57 + t58 * t65 + t99, t54 * t89 + t61 * t65 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t17 + t96, qJD(3) * t18 - t97, t53, t43, 0, 0, 0, qJD(3) * t3 + t36 * t92, -qJD(3) * t2 - t77; 0, 0, 0, 0, 0, 0, 0, 0, t83 + t87, t82 + t86, t53, t43, 0, 0, 0, t61 * t83 + t100, -pkin(3) * t75 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t43, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t115, t57, -t89, 0, -pkin(6) * t57, pkin(6) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t115, 0, 0, 0, -qJD(3) * t5 - t79, -qJD(3) * t7 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t115, 0, 0, 0, t58 * t82 - t99, t61 * t82 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t115, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
