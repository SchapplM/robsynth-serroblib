% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:53
% DurationCPUTime: 0.44s
% Computational Cost: add. (1396->97), mult. (1997->155), div. (0->0), fcn. (1548->10), ass. (0->77)
t70 = qJD(3) + qJD(4);
t68 = t70 ^ 2;
t69 = qJDD(3) + qJDD(4);
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t49 = t68 * t79 - t69 * t82;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t87 = -t68 * t82 - t69 * t79;
t99 = t49 * t80 + t83 * t87;
t98 = t49 * t83 - t80 * t87;
t75 = sin(pkin(8));
t77 = cos(pkin(8));
t58 = -g(1) * t77 - g(2) * t75;
t73 = -g(3) + qJDD(1);
t74 = sin(pkin(9));
t76 = cos(pkin(9));
t38 = -t58 * t74 + t73 * t76;
t39 = t58 * t76 + t73 * t74;
t25 = t38 * t80 + t39 * t83;
t85 = qJD(3) ^ 2;
t19 = -pkin(3) * t85 + t25;
t24 = t83 * t38 - t39 * t80;
t86 = qJDD(3) * pkin(3) + t24;
t17 = t82 * t19 + t79 * t86;
t15 = -pkin(4) * t68 + pkin(7) * t69 + t17;
t55 = -g(1) * t75 + g(2) * t77 + qJDD(2);
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t10 = t15 * t81 + t55 * t78;
t9 = t15 * t78 - t55 * t81;
t4 = t81 * t10 + t78 * t9;
t16 = -t19 * t79 + t82 * t86;
t14 = -pkin(4) * t69 - pkin(7) * t68 - t16;
t95 = -pkin(4) * t14 + pkin(7) * t4;
t61 = t78 * t68 * t81;
t53 = qJDD(5) + t61;
t94 = t78 * t53;
t93 = t78 * t69;
t54 = qJDD(5) - t61;
t92 = t81 * t54;
t91 = qJD(5) * t70;
t71 = t78 ^ 2;
t64 = t71 * t68;
t84 = qJD(5) ^ 2;
t59 = -t64 - t84;
t36 = -t59 * t78 - t92;
t43 = 0.2e1 * t81 * t91 + t93;
t90 = -pkin(4) * t43 + pkin(7) * t36 + t78 * t14;
t72 = t81 ^ 2;
t65 = t72 * t68;
t60 = -t65 - t84;
t35 = t60 * t81 - t94;
t63 = t81 * t69;
t44 = -0.2e1 * t78 * t91 + t63;
t89 = pkin(4) * t44 + pkin(7) * t35 - t81 * t14;
t46 = (t71 + t72) * t69;
t52 = t64 + t65;
t88 = pkin(4) * t52 + pkin(7) * t46 + t4;
t57 = -qJDD(3) * t80 - t83 * t85;
t56 = qJDD(3) * t83 - t80 * t85;
t34 = t94 + t81 * (-t64 + t84);
t33 = t78 * (t65 - t84) + t92;
t30 = t43 * t78;
t29 = t44 * t81;
t28 = t46 * t82 - t52 * t79;
t27 = t46 * t79 + t52 * t82;
t26 = t43 * t81 + t44 * t78;
t23 = t36 * t82 + t43 * t79;
t22 = t35 * t82 - t44 * t79;
t21 = t36 * t79 - t43 * t82;
t20 = t35 * t79 + t44 * t82;
t6 = -t16 * t79 + t17 * t82;
t5 = t16 * t82 + t17 * t79;
t2 = t14 * t79 + t4 * t82;
t1 = -t14 * t82 + t4 * t79;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t76 + t39 * t74, 0, 0, 0, 0, 0, 0, t56 * t76 + t57 * t74, -t56 * t74 + t57 * t76, 0, t74 * (-t24 * t80 + t25 * t83) + t76 * (t24 * t83 + t25 * t80), 0, 0, 0, 0, 0, 0, t74 * t99 - t76 * t98, t74 * t98 + t76 * t99, 0, t74 * (-t5 * t80 + t6 * t83) + t76 * (t5 * t83 + t6 * t80), 0, 0, 0, 0, 0, 0, t74 * (-t20 * t80 + t22 * t83) + t76 * (t20 * t83 + t22 * t80), t74 * (-t21 * t80 + t23 * t83) + t76 * (t21 * t83 + t23 * t80), t74 * (-t27 * t80 + t28 * t83) + t76 * (t27 * t83 + t28 * t80), t74 * (-t1 * t80 + t2 * t83) + t76 * (t1 * t83 + t2 * t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, t53 * t81 + t60 * t78, -t54 * t78 + t59 * t81, 0, t10 * t78 - t81 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t24, -t25, 0, 0, 0, 0, 0, 0, 0, t69, -pkin(3) * t49 + t16, pkin(3) * t87 - t17, 0, pkin(3) * t5, t30, t26, t34, t29, t33, 0, pkin(3) * t20 + t89, pkin(3) * t21 + t90, pkin(3) * t27 + t88, pkin(3) * t1 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t16, -t17, 0, 0, t30, t26, t34, t29, t33, 0, t89, t90, t88, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t64 - t65, t93, t61, t63, qJDD(5), -t9, -t10, 0, 0;];
tauJ_reg = t3;
