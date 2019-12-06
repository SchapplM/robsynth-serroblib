% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:37
% EndTime: 2019-12-05 15:47:40
% DurationCPUTime: 0.49s
% Computational Cost: add. (506->70), mult. (1137->111), div. (0->0), fcn. (1262->8), ass. (0->67)
t74 = qJD(4) + qJD(5);
t64 = cos(qJ(5));
t65 = cos(qJ(4));
t96 = t64 * t65;
t62 = sin(qJ(5));
t63 = sin(qJ(4));
t99 = t62 * t63;
t51 = -t96 + t99;
t106 = t74 * t51;
t97 = t64 * t63;
t98 = t62 * t65;
t53 = t97 + t98;
t31 = t74 * t53;
t91 = sin(pkin(9));
t59 = t91 * pkin(2) + pkin(6);
t102 = pkin(7) + t59;
t47 = t102 * t63;
t48 = t102 * t65;
t100 = sin(qJ(2));
t101 = cos(qJ(2));
t92 = cos(pkin(9));
t49 = t91 * t100 - t92 * t101;
t67 = t98 / 0.2e1 + t97 / 0.2e1;
t9 = (-t53 / 0.2e1 + t67) * t49;
t93 = t9 * qJD(1);
t105 = -t93 + t74 * (t62 * t47 - t64 * t48);
t66 = t96 / 0.2e1 - t99 / 0.2e1;
t10 = (t51 / 0.2e1 + t66) * t49;
t85 = t10 * qJD(1);
t104 = -t85 + t74 * (t64 * t47 + t62 * t48);
t103 = pkin(4) * t63;
t95 = pkin(4) * qJD(5);
t94 = qJD(4) * pkin(4);
t90 = qJD(2) * t53;
t60 = -t92 * pkin(2) - pkin(3);
t54 = -t65 * pkin(4) + t60;
t89 = qJD(2) * t54;
t88 = qJD(2) * t63;
t87 = qJD(2) * t65;
t86 = qJD(5) * t54;
t18 = t51 * t103 + t54 * t53;
t80 = t18 * qJD(2);
t19 = t53 * t103 - t54 * t51;
t79 = t19 * qJD(2);
t20 = t51 ^ 2 - t53 ^ 2;
t78 = t20 * qJD(2);
t55 = -t63 ^ 2 + t65 ^ 2;
t77 = t55 * qJD(2);
t76 = t63 * qJD(4);
t75 = t65 * qJD(4);
t73 = t51 * t89;
t72 = t53 * t89;
t71 = t60 * t88;
t70 = t60 * t87;
t69 = t63 * t87;
t68 = pkin(4) * t74;
t50 = -t92 * t100 - t91 * t101;
t34 = t51 * t90;
t28 = t49 * t65;
t27 = t49 * t63;
t12 = (t53 / 0.2e1 + t67) * t49;
t11 = (-t51 / 0.2e1 + t66) * t49;
t5 = t10 * qJD(2);
t3 = t9 * qJD(2);
t2 = t12 * qJD(2) - t106 * t50;
t1 = t11 * qJD(2) - t50 * t31;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t100 * qJD(2), -t101 * qJD(2), (-t91 * t49 + t92 * t50) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, t27 * qJD(4) + t50 * t87, t28 * qJD(4) - t50 * t88, 0, 0, 0, 0, 0, -t50 * t51 * qJD(2) + t74 * t12, t74 * t11 - t50 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * qJD(2) + t50 * t75, t28 * qJD(2) - t50 * t76, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t9, -t74 * t10; 0, 0, 0, 0, 0, t63 * t75, t55 * qJD(4), 0, 0, 0, t60 * t76, t60 * t75, -t51 * t31, t74 * t20, 0, 0, 0, t18 * qJD(4) + t53 * t86, t19 * qJD(4) - t51 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t69, t77, t75, -t76, 0, -t59 * t75 + t71, t59 * t76 + t70, -t34, t78, -t106, -t31, 0, t80 + t105, t79 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t78, -t106, -t31, 0, t72 + t105, -t73 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0, 0, 0, 0, -t31, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5; 0, 0, 0, 0, 0, -t69, -t77, 0, 0, 0, -t71, -t70, t34, -t78, 0, 0, 0, t93 - t80, t85 - t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t95, -t64 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t68, -t64 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t78, 0, 0, 0, t93 - t72, t85 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t94, t64 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
