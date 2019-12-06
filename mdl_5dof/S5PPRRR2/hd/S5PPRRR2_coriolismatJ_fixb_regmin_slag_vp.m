% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:52
% EndTime: 2019-12-05 15:14:54
% DurationCPUTime: 0.43s
% Computational Cost: add. (433->66), mult. (1043->104), div. (0->0), fcn. (1154->8), ass. (0->65)
t76 = qJD(4) + qJD(5);
t58 = sin(qJ(4));
t56 = sin(pkin(9));
t59 = sin(qJ(3));
t92 = cos(pkin(9));
t96 = cos(qJ(3));
t70 = t59 * t56 - t96 * t92;
t68 = t70 * t58;
t65 = t68 / 0.2e1;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t61 = cos(qJ(4));
t47 = t57 * t58 - t60 * t61;
t101 = t76 * t47;
t49 = t57 * t61 + t60 * t58;
t28 = t76 * t49;
t98 = pkin(6) + pkin(7);
t50 = t98 * t58;
t51 = t98 * t61;
t64 = t70 * t61 / 0.2e1;
t62 = t57 * t64 + t60 * t65;
t66 = t70 * t49;
t9 = -t66 / 0.2e1 + t62;
t93 = t9 * qJD(1);
t100 = -t93 + t76 * (t57 * t50 - t60 * t51);
t63 = -t57 * t68 / 0.2e1 + t60 * t64;
t67 = t70 * t47;
t10 = t67 / 0.2e1 + t63;
t88 = t10 * qJD(1);
t99 = -t88 + t76 * (t60 * t50 + t57 * t51);
t97 = pkin(4) * t58;
t95 = pkin(4) * qJD(5);
t94 = qJD(4) * pkin(4);
t55 = -t61 * pkin(4) - pkin(3);
t91 = qJD(3) * t55;
t90 = qJD(3) * t61;
t89 = qJD(5) * t55;
t15 = t47 ^ 2 - t49 ^ 2;
t87 = t15 * qJD(3);
t20 = t47 * t97 + t55 * t49;
t82 = t20 * qJD(3);
t21 = -t55 * t47 + t49 * t97;
t81 = t21 * qJD(3);
t46 = t96 * t56 + t59 * t92;
t80 = t46 * qJD(3);
t52 = -t58 ^ 2 + t61 ^ 2;
t79 = t52 * qJD(3);
t78 = t58 * qJD(4);
t77 = t61 * qJD(4);
t75 = pkin(3) * t58 * qJD(3);
t74 = pkin(3) * t90;
t73 = t47 * t91;
t72 = t49 * t91;
t71 = t58 * t90;
t69 = pkin(4) * t76;
t31 = t49 * t47 * qJD(3);
t25 = 0.2e1 * t64;
t24 = 0.2e1 * t65;
t12 = t66 / 0.2e1 + t62;
t11 = -t67 / 0.2e1 + t63;
t5 = t10 * qJD(3);
t3 = t9 * qJD(3);
t2 = t12 * qJD(3) + t101 * t46;
t1 = t11 * qJD(3) + t46 * t28;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t80, t70 * qJD(3), 0, 0, 0, 0, 0, t24 * qJD(4) - t61 * t80, t25 * qJD(4) + t58 * t80, 0, 0, 0, 0, 0, t76 * t12 + t47 * t80, t76 * t11 + t49 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * qJD(3) - t46 * t77, t25 * qJD(3) + t46 * t78, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77, 0, 0, 0, 0, 0, -t28, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76 * t9, -t76 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t58 * t77, t52 * qJD(4), 0, 0, 0, -pkin(3) * t78, -pkin(3) * t77, -t47 * t28, t76 * t15, 0, 0, 0, t20 * qJD(4) + t49 * t89, t21 * qJD(4) - t47 * t89; 0, 0, 0, 0, 0, t71, t79, t77, -t78, 0, -pkin(6) * t77 - t75, pkin(6) * t78 - t74, -t31, t87, -t101, -t28, 0, t82 + t100, t81 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t87, -t101, -t28, 0, t72 + t100, -t73 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t71, -t79, 0, 0, 0, t75, t74, t31, -t87, 0, 0, 0, t93 - t82, t88 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t95, -t60 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t69, -t60 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t87, 0, 0, 0, t93 - t72, t88 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t94, t60 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
