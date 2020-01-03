% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:43
% DurationCPUTime: 0.49s
% Computational Cost: add. (850->117), mult. (1934->180), div. (0->0), fcn. (1326->8), ass. (0->83)
t69 = sin(pkin(6));
t71 = cos(pkin(6));
t55 = -t71 * g(1) - t69 * g(2);
t73 = sin(qJ(2));
t75 = cos(qJ(2));
t90 = -g(3) + qJDD(1);
t37 = t75 * t55 + t73 * t90;
t77 = qJD(2) ^ 2;
t83 = -t77 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t37;
t68 = sin(pkin(7));
t70 = cos(pkin(7));
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t46 = (t68 * t72 - t70 * t74) * qJD(2);
t97 = t70 * t72;
t82 = t68 * t74 + t97;
t48 = t82 * qJD(2);
t32 = t48 * t46;
t101 = qJDD(4) - t32;
t105 = t101 * t72;
t104 = t101 * t74;
t103 = t70 * t77;
t65 = t70 ^ 2;
t78 = t68 ^ 2;
t102 = t65 + t78;
t53 = t102 * t77;
t43 = t46 ^ 2;
t44 = t48 ^ 2;
t54 = -t69 * g(1) + t71 * g(2);
t86 = t68 * t54 + t83 * t70;
t88 = t70 * qJDD(2);
t99 = t65 * t77;
t16 = -pkin(3) * t99 + pkin(5) * t88 + t86;
t81 = (pkin(3) * t103 - pkin(5) * qJDD(2) - t83) * t68;
t98 = t70 * t54;
t7 = t72 * t16 - t74 * (t81 + t98);
t8 = t74 * t16 + t54 * t97 + t72 * t81;
t2 = -t74 * t7 + t72 * t8;
t100 = t68 * t2;
t36 = -t73 * t55 + t75 * t90;
t67 = qJDD(2) * pkin(2);
t34 = -t77 * qJ(3) + qJDD(3) - t36 - t67;
t23 = -pkin(3) * t88 + t34 + (-t78 * t77 - t99) * pkin(5);
t96 = t72 * t23;
t26 = qJDD(4) + t32;
t95 = t72 * t26;
t94 = t74 * t23;
t93 = t74 * t26;
t92 = t46 * qJD(4);
t91 = t48 * qJD(4);
t89 = t68 * qJDD(2);
t87 = t75 * qJDD(2);
t3 = t72 * t7 + t74 * t8;
t10 = t68 * (t83 * t68 - t98) + t70 * t86;
t84 = -t34 + t67;
t20 = -t72 * t89 + t74 * t88;
t45 = t82 * qJDD(2);
t76 = qJD(4) ^ 2;
t64 = t65 * qJDD(2);
t63 = t78 * qJDD(2);
t52 = t64 + t63;
t50 = t102 * t103;
t49 = t68 * t53;
t40 = -t44 - t76;
t39 = -t44 + t76;
t38 = t43 - t76;
t31 = t45 - t92;
t30 = t45 - 0.2e1 * t92;
t29 = t20 - t91;
t28 = -t20 + 0.2e1 * t91;
t24 = -t76 - t43;
t22 = -t43 - t44;
t18 = -t72 * t40 - t93;
t17 = t74 * t40 - t95;
t14 = t74 * t20 + t72 * t45;
t13 = t72 * t20 - t74 * t45;
t12 = t74 * t24 - t105;
t11 = t72 * t24 + t104;
t9 = -t68 * t17 + t70 * t18;
t5 = -t68 * t13 + t70 * t14;
t4 = -t68 * t11 + t70 * t12;
t1 = t70 * t3 - t100;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, 0, 0, 0, 0, 0, -t73 * t77 + t87, -t73 * qJDD(2) - t75 * t77, 0, t75 * t36 + t73 * t37, 0, 0, 0, 0, 0, 0, -t73 * t50 + t70 * t87, t73 * t49 - t68 * t87, t73 * t52 + t75 * t53, t73 * t10 - t75 * t34, 0, 0, 0, 0, 0, 0, -t75 * t28 + t73 * t4, -t75 * t30 + t73 * t9, -t75 * t22 + t73 * t5, t73 * t1 - t75 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t36, -t37, 0, 0, t63, 0.2e1 * t68 * t88, 0, t64, 0, 0, -qJ(3) * t50 + t84 * t70, qJ(3) * t49 - t84 * t68, pkin(2) * t53 + qJ(3) * t52 + t10, -pkin(2) * t34 + qJ(3) * t10, t68 * (t74 * t31 - t72 * t91) + t70 * (t72 * t31 + t74 * t91), t68 * (-t74 * t28 - t72 * t30) + t70 * (-t72 * t28 + t74 * t30), t68 * (-t72 * t39 + t104) + t70 * (t74 * t39 + t105), t68 * (-t72 * t29 + t74 * t92) + t70 * (t74 * t29 + t72 * t92), t68 * (t74 * t38 - t95) + t70 * (t72 * t38 + t93), (t68 * (-t46 * t74 + t48 * t72) + t70 * (-t46 * t72 - t48 * t74)) * qJD(4), t68 * (-pkin(5) * t11 + t96) + t70 * (-pkin(3) * t28 + pkin(5) * t12 - t94) - pkin(2) * t28 + qJ(3) * t4, t68 * (-pkin(5) * t17 + t94) + t70 * (-pkin(3) * t30 + pkin(5) * t18 + t96) - pkin(2) * t30 + qJ(3) * t9, t68 * (-pkin(5) * t13 - t2) + t70 * (-pkin(3) * t22 + pkin(5) * t14 + t3) - pkin(2) * t22 + qJ(3) * t5, -pkin(5) * t100 + t70 * (-pkin(3) * t23 + pkin(5) * t3) - pkin(2) * t23 + qJ(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t89, -t53, t34, 0, 0, 0, 0, 0, 0, t28, t30, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t44 - t43, t45, -t32, t20, qJDD(4), -t7, -t8, 0, 0;];
tauJ_reg = t6;
