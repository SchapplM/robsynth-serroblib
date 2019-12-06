% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:07
% EndTime: 2019-12-05 15:31:10
% DurationCPUTime: 0.45s
% Computational Cost: add. (557->76), mult. (1249->133), div. (0->0), fcn. (1058->4), ass. (0->81)
t52 = sin(qJ(4));
t48 = t52 ^ 2;
t53 = cos(qJ(4));
t49 = t53 ^ 2;
t101 = t48 + t49;
t99 = t52 * pkin(4);
t51 = cos(pkin(8));
t50 = sin(pkin(8));
t34 = -pkin(3) * t51 - pkin(6) * t50 - pkin(2);
t27 = t53 * t34;
t94 = qJ(5) * t50;
t72 = t53 * t94;
t61 = t27 - t72;
t96 = qJ(3) * t52;
t18 = (-pkin(4) - t96) * t51 + t61;
t98 = t18 * t52;
t97 = t50 * t53;
t46 = t50 ^ 2;
t39 = t51 ^ 2 + t46;
t95 = qJ(3) * t53;
t56 = -t27 / 0.2e1 + t18 / 0.2e1 + t72 / 0.2e1;
t66 = t96 / 0.2e1;
t2 = ((t66 + pkin(4) / 0.2e1) * t51 + t56) * t97;
t93 = t2 * qJD(2);
t73 = t51 * t96;
t55 = t61 - t73;
t54 = t55 * t52;
t3 = (t54 - t98) * t50;
t92 = t3 * qJD(2);
t4 = ((-pkin(4) / 0.2e1 + t66) * t51 + t56) * t52;
t91 = t4 * qJD(2);
t24 = -t52 * t34 - t51 * t95;
t21 = -t52 * t94 - t24;
t8 = (t18 * t53 + t21 * t52) * t50;
t90 = t8 * qJD(2);
t20 = t24 * t51 - t46 * t95;
t89 = qJD(2) * t20;
t88 = qJD(2) * t51;
t87 = qJD(3) * t51;
t86 = qJD(4) * t52;
t85 = qJD(4) * t53;
t23 = -t27 + t73;
t19 = -t23 * t51 - t46 * t96;
t84 = t19 * qJD(2);
t67 = t48 / 0.2e1 + t49 / 0.2e1;
t22 = (-0.1e1 / 0.2e1 + t67) * t51 * t50;
t83 = t22 * qJD(2);
t25 = (0.1e1 / 0.2e1 + t67) * t50;
t82 = t25 * qJD(2);
t28 = t101 * t46;
t81 = t28 * qJD(2);
t29 = (t48 - t49) * t46;
t80 = t29 * qJD(2);
t30 = t39 * t52;
t79 = t30 * qJD(2);
t31 = t39 * t53;
t78 = t31 * qJD(2);
t37 = t39 * qJ(3);
t77 = t37 * qJD(2);
t76 = t39 * qJD(2);
t75 = pkin(4) * t97;
t74 = t46 * t52 * t53;
t71 = t50 * t86;
t70 = t50 * t85;
t69 = t52 * t88;
t68 = t53 * t88;
t65 = pkin(4) * t70;
t64 = qJD(2) * t75;
t63 = qJD(2) * t74;
t62 = -qJD(4) + t88;
t32 = (qJ(3) + t99) * t50;
t6 = t32 * t75 + (-t18 + t55) * t21;
t60 = -t2 * qJD(1) + t6 * qJD(2);
t7 = t32 * t50 + (t21 * t53 - t98) * t51;
t59 = -qJD(1) * t22 - qJD(2) * t7;
t58 = t62 * t52;
t57 = t62 * t53;
t26 = (0.1e1 / 0.2e1 - t101 / 0.2e1) * t50;
t5 = t54 / 0.2e1 - t98 / 0.2e1 - t51 * t99 / 0.2e1;
t1 = qJD(3) * t22 - qJD(4) * t2;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t71, 0, -t65 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, t39 * qJD(3), t37 * qJD(3), -qJD(4) * t74, t29 * qJD(4), t51 * t71, t51 * t70, 0, qJD(3) * t30 - qJD(4) * t20, qJD(3) * t31 + qJD(4) * t19, -qJD(4) * t3 + qJD(5) * t28, qJD(3) * t7 + qJD(4) * t6 - qJD(5) * t8; 0, 0, 0, 0, 0, 0, t76, t77, 0, 0, 0, 0, 0, t79, t78, 0, qJD(4) * t5 + qJD(5) * t26 - t59; 0, 0, 0, 0, 0, 0, 0, 0, -t63, t80, t50 * t58, t50 * t57, 0, qJD(4) * t24 - t89, qJD(4) * t23 + t84, pkin(4) * t71 - t92, -pkin(4) * qJD(4) * t21 + qJD(3) * t5 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, qJD(3) * t26 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83; 0, 0, 0, 0, 0, 0, -t76, -t77, 0, 0, 0, 0, 0, t51 * t86 - t79, t51 * t85 - t78, 0, -qJD(4) * t4 - qJD(5) * t25 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57, 0, -pkin(4) * t86 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, t63, -t80, -t50 * t69, -t50 * t68, 0, -t52 * t87 + t89, -t53 * t87 - t84, t92, t4 * qJD(3) - qJD(5) * t75 - t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68, 0, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t25 * qJD(3) + t65 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
