% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:19
% EndTime: 2019-12-05 17:06:21
% DurationCPUTime: 0.46s
% Computational Cost: add. (2031->92), mult. (2834->134), div. (0->0), fcn. (2000->10), ass. (0->79)
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t60 = g(1) * t77 - g(2) * t78;
t61 = -g(1) * t78 - g(2) * t77;
t82 = sin(qJ(2));
t86 = cos(qJ(2));
t95 = t86 * t60 - t61 * t82;
t39 = qJDD(2) * pkin(2) + t95;
t90 = -t60 * t82 - t61 * t86;
t40 = -qJD(2) ^ 2 * pkin(2) - t90;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t21 = t85 * t39 - t81 * t40;
t72 = qJDD(2) + qJDD(3);
t19 = pkin(3) * t72 + t21;
t22 = t81 * t39 + t85 * t40;
t73 = qJD(2) + qJD(3);
t71 = t73 ^ 2;
t20 = -pkin(3) * t71 + t22;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t17 = t80 * t19 + t84 * t20;
t70 = qJD(4) + t73;
t68 = t70 ^ 2;
t69 = qJDD(4) + t72;
t15 = -pkin(4) * t68 + pkin(8) * t69 + t17;
t76 = -g(3) + qJDD(1);
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t10 = t15 * t83 + t76 * t79;
t9 = t15 * t79 - t76 * t83;
t4 = t83 * t10 + t79 * t9;
t16 = t84 * t19 - t80 * t20;
t14 = -pkin(4) * t69 - pkin(8) * t68 - t16;
t104 = -pkin(4) * t14 + pkin(8) * t4;
t58 = t79 * t68 * t83;
t53 = qJDD(5) + t58;
t103 = t79 * t53;
t102 = t79 * t69;
t54 = qJDD(5) - t58;
t101 = t83 * t54;
t100 = qJD(5) * t70;
t2 = -t14 * t84 + t4 * t80;
t99 = pkin(3) * t2 + t104;
t74 = t79 ^ 2;
t63 = t74 * t68;
t87 = qJD(5) ^ 2;
t56 = -t63 - t87;
t38 = -t56 * t79 - t101;
t46 = 0.2e1 * t100 * t83 + t102;
t98 = -pkin(4) * t46 + pkin(8) * t38 + t79 * t14;
t75 = t83 ^ 2;
t64 = t75 * t68;
t57 = -t64 - t87;
t36 = t57 * t83 - t103;
t62 = t83 * t69;
t47 = -0.2e1 * t100 * t79 + t62;
t97 = pkin(4) * t47 + pkin(8) * t36 - t83 * t14;
t49 = (t74 + t75) * t69;
t52 = t63 + t64;
t96 = pkin(4) * t52 + pkin(8) * t49 + t4;
t26 = t38 * t80 - t46 * t84;
t94 = pkin(3) * t26 + t98;
t25 = t36 * t80 + t47 * t84;
t93 = pkin(3) * t25 + t97;
t29 = t49 * t80 + t52 * t84;
t92 = pkin(3) * t29 + t96;
t89 = t68 * t80 - t69 * t84;
t91 = -pkin(3) * t89 + t16;
t50 = -t68 * t84 - t69 * t80;
t88 = pkin(3) * t50 - t17;
t37 = t101 + t79 * (t64 - t87);
t35 = t83 * (-t63 + t87) + t103;
t31 = t47 * t83;
t30 = t46 * t79;
t27 = t46 * t83 + t47 * t79;
t6 = t16 * t84 + t17 * t80;
t5 = pkin(3) * t6;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, t53 * t83 + t57 * t79, -t54 * t79 + t56 * t83, 0, t10 * t79 - t83 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t95, t90, 0, 0, 0, 0, 0, 0, 0, t72, pkin(2) * (-t71 * t81 + t72 * t85) + t21, pkin(2) * (-t71 * t85 - t72 * t81) - t22, 0, pkin(2) * (t21 * t85 + t22 * t81), 0, 0, 0, 0, 0, t69, pkin(2) * (t50 * t81 - t85 * t89) + t91, pkin(2) * (t85 * t50 + t81 * t89) + t88, 0, pkin(2) * (t81 * (-t16 * t80 + t17 * t84) + t85 * t6) + t5, t30, t27, t35, t31, t37, 0, pkin(2) * (t81 * (t36 * t84 - t47 * t80) + t85 * t25) + t93, pkin(2) * (t81 * (t38 * t84 + t46 * t80) + t85 * t26) + t94, pkin(2) * (t81 * (t49 * t84 - t52 * t80) + t85 * t29) + t92, pkin(2) * (t81 * (t14 * t80 + t4 * t84) + t85 * t2) + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t21, -t22, 0, 0, 0, 0, 0, 0, 0, t69, t91, t88, 0, t5, t30, t27, t35, t31, t37, 0, t93, t94, t92, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t16, -t17, 0, 0, t30, t27, t35, t31, t37, 0, t97, t98, t96, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t63 - t64, t102, t58, t62, qJDD(5), -t9, -t10, 0, 0;];
tauJ_reg = t1;
