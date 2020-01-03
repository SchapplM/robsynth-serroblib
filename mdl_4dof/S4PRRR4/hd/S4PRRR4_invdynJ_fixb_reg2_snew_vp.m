% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRR4
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:41
% DurationCPUTime: 0.56s
% Computational Cost: add. (1220->144), mult. (2531->195), div. (0->0), fcn. (1700->8), ass. (0->95)
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t77 = cos(qJ(3));
t74 = sin(qJ(3));
t95 = qJD(2) * t74;
t41 = -t76 * t77 * qJD(2) + t73 * t95;
t43 = (t77 * t73 + t74 * t76) * qJD(2);
t34 = t43 * t41;
t66 = qJDD(3) + qJDD(4);
t110 = -t34 + t66;
t112 = t110 * t73;
t111 = t110 * t76;
t79 = qJD(2) ^ 2;
t107 = cos(qJ(2));
t96 = sin(pkin(7));
t97 = cos(pkin(7));
t51 = -t97 * g(1) - t96 * g(2);
t75 = sin(qJ(2));
t82 = t96 * g(1) - t97 * g(2);
t81 = -t107 * t51 - t75 * t82;
t94 = qJDD(2) * pkin(5);
t32 = -t79 * pkin(2) - t81 + t94;
t70 = -g(3) + qJDD(1);
t25 = t74 * t32 - t77 * t70;
t61 = t74 * qJDD(2);
t91 = qJD(2) * qJD(3);
t89 = t77 * t91;
t48 = t61 + t89;
t109 = -t25 + (-t48 + t89) * pkin(6);
t39 = t41 ^ 2;
t40 = t43 ^ 2;
t67 = qJD(3) + qJD(4);
t65 = t67 ^ 2;
t26 = t77 * t32 + t74 * t70;
t62 = t77 * qJDD(2);
t90 = t74 * t91;
t49 = t62 - t90;
t54 = qJD(3) * pkin(3) - pkin(6) * t95;
t69 = t77 ^ 2;
t64 = t69 * t79;
t11 = -pkin(3) * t64 + t49 * pkin(6) - qJD(3) * t54 + t26;
t57 = t74 * t79 * t77;
t92 = qJDD(3) + t57;
t80 = t92 * pkin(3) + t109;
t4 = t73 * t11 - t76 * t80;
t101 = t76 * t11;
t5 = t73 * t80 + t101;
t1 = -t76 * t4 + t73 * t5;
t108 = t74 * t1;
t106 = t67 * t73;
t105 = t67 * t76;
t86 = t107 * t82 - t75 * t51;
t31 = -qJDD(2) * pkin(2) - t79 * pkin(5) - t86;
t20 = -t49 * pkin(3) - pkin(6) * t64 + t54 * t95 + t31;
t104 = t73 * t20;
t29 = t34 + t66;
t103 = t73 * t29;
t102 = t74 * t92;
t100 = t76 * t20;
t99 = t76 * t29;
t53 = qJDD(3) - t57;
t98 = t77 * t53;
t93 = qJD(4) + t67;
t2 = t73 * t4 + t76 * t5;
t88 = t74 * t25 + t77 * t26;
t87 = t73 * t48 - t76 * t49;
t84 = t76 * t48 + t73 * t49;
t83 = (-qJD(4) + t67) * t43 - t87;
t22 = -t41 * qJD(4) + t84;
t78 = qJD(3) ^ 2;
t68 = t74 ^ 2;
t63 = t68 * t79;
t56 = -t64 - t78;
t55 = -t63 - t78;
t50 = t62 - 0.2e1 * t90;
t47 = t61 + 0.2e1 * t89;
t38 = t67 * t41;
t37 = -t40 + t65;
t36 = t39 - t65;
t35 = -t40 - t65;
t33 = t40 - t39;
t27 = -t65 - t39;
t23 = -t39 - t40;
t21 = -t43 * qJD(4) - t87;
t19 = -t73 * t35 - t99;
t18 = t76 * t35 - t103;
t17 = t22 + t38;
t16 = t22 - t38;
t15 = -t93 * t41 + t84;
t12 = t93 * t43 + t87;
t9 = t76 * t27 - t112;
t8 = t73 * t27 + t111;
t7 = t73 * t17 + t76 * t83;
t6 = -t76 * t17 + t73 * t83;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, t74 * t56 + t77 * t92, -t74 * t53 + t77 * t55, 0, -t77 * t25 + t74 * t26, 0, 0, 0, 0, 0, 0, t74 * t9 + t77 * t8, t77 * t18 + t74 * t19, t77 * t6 + t74 * t7, t77 * t1 + t74 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t86, t81, 0, 0, (t48 + t89) * t74, t77 * t47 + t74 * t50, t102 + t77 * (-t63 + t78), (t49 - t90) * t77, t74 * (t64 - t78) + t98, 0, -t77 * t31 + pkin(2) * t50 + pkin(5) * (t77 * t56 - t102), t74 * t31 - pkin(2) * t47 + pkin(5) * (-t74 * t55 - t98), pkin(2) * (t63 + t64) + (t68 + t69) * t94 + t88, -pkin(2) * t31 + pkin(5) * t88, t74 * (-t43 * t106 + t76 * t22) + t77 * (t43 * t105 + t73 * t22), t74 * (-t76 * t12 - t73 * t16) + t77 * (-t73 * t12 + t76 * t16), t74 * (-t73 * t37 + t111) + t77 * (t76 * t37 + t112), t74 * (t41 * t105 - t73 * t21) + t77 * (t41 * t106 + t76 * t21), t74 * (t76 * t36 - t103) + t77 * (t73 * t36 + t99), (t74 * (-t41 * t76 + t43 * t73) + t77 * (-t41 * t73 - t43 * t76)) * t67, t74 * (-pkin(6) * t8 + t104) + t77 * (-pkin(3) * t12 + pkin(6) * t9 - t100) - pkin(2) * t12 + pkin(5) * (-t74 * t8 + t77 * t9), t74 * (-pkin(6) * t18 + t100) + t77 * (-pkin(3) * t15 + pkin(6) * t19 + t104) - pkin(2) * t15 + pkin(5) * (-t74 * t18 + t77 * t19), t74 * (-pkin(6) * t6 - t1) + t77 * (-pkin(3) * t23 + pkin(6) * t7 + t2) - pkin(2) * t23 + pkin(5) * (-t74 * t6 + t77 * t7), -pkin(6) * t108 + t77 * (-pkin(3) * t20 + pkin(6) * t2) - pkin(2) * t20 + pkin(5) * (t77 * t2 - t108); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t63 - t64, t61, t57, t62, qJDD(3), -t25, -t26, 0, 0, t34, t33, t17, -t34, t83, t66, pkin(3) * t8 - t4, -t101 - t73 * t109 + (-t73 * t92 + t18) * pkin(3), pkin(3) * t6, pkin(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, t17, -t34, t83, t66, -t4, -t5, 0, 0;];
tauJ_reg = t3;
