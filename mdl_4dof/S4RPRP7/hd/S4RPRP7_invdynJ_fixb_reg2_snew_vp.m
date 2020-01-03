% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:18
% DurationCPUTime: 0.37s
% Computational Cost: add. (435->104), mult. (892->102), div. (0->0), fcn. (396->4), ass. (0->79)
t52 = cos(qJ(3));
t50 = sin(qJ(3));
t67 = pkin(3) * t52 + qJ(4) * t50;
t99 = t67 * qJD(3) - 0.2e1 * qJD(4) * t52;
t94 = -pkin(5) - pkin(1);
t54 = qJD(3) ^ 2;
t49 = t52 ^ 2;
t55 = qJD(1) ^ 2;
t89 = t49 * t55;
t37 = t54 + t89;
t38 = t50 * t55 * t52;
t33 = qJDD(3) + t38;
t88 = t50 * t33;
t11 = t37 * t52 + t88;
t98 = t94 * t11;
t77 = qJD(1) * qJD(3);
t72 = t52 * t77;
t79 = t50 * qJDD(1);
t27 = -t72 - t79;
t73 = t50 * t77;
t78 = t52 * qJDD(1);
t28 = -t73 + t78;
t97 = t27 * pkin(3) + t28 * qJ(4);
t47 = qJDD(1) * qJ(2);
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t85 = t53 * g(1) + t51 * g(2);
t75 = -t47 + t85;
t61 = -t94 * t55 + t75;
t95 = t97 + (-(2 * qJD(2)) - t99) * qJD(1);
t93 = t50 * g(3);
t92 = t52 * g(3);
t91 = t54 * pkin(3);
t48 = t50 ^ 2;
t90 = t48 * t55;
t74 = t51 * g(1) - t53 * g(2);
t69 = qJDD(2) - t74;
t59 = -t55 * qJ(2) + t69;
t17 = qJDD(1) * t94 + t59;
t87 = t52 * t17;
t34 = qJDD(3) - t38;
t86 = t52 * t34;
t84 = t48 + t49;
t83 = t52 * qJ(4);
t66 = pkin(3) * t50 - t83;
t82 = t55 * t66;
t80 = qJDD(1) * pkin(1);
t76 = qJD(2) * qJD(1);
t36 = -t54 - t90;
t10 = t36 * t50 + t86;
t26 = 0.2e1 * t72 + t79;
t71 = qJ(2) * t26 + t10 * t94;
t30 = t84 * qJDD(1);
t32 = t84 * t55;
t70 = -qJ(2) * t32 - t30 * t94;
t6 = t87 + t93;
t15 = t50 * t17;
t7 = -t15 + t92;
t3 = -t50 * t7 + t52 * t6;
t68 = -t52 * t82 - qJDD(4) + t87;
t29 = -0.2e1 * t73 + t78;
t65 = t26 * t52 + t29 * t50;
t64 = -t52 * (-t54 + t90) + t88;
t63 = -t55 * pkin(1) - t75;
t62 = qJ(2) + t66;
t60 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t50 * t82 + t15;
t58 = t60 - t92;
t57 = -qJDD(3) * pkin(3) - t54 * qJ(4) - t68;
t42 = 0.2e1 * t76;
t31 = (-t48 + t49) * t55;
t18 = -t59 + t80;
t16 = t61 - 0.2e1 * t76;
t14 = t86 - t50 * (t54 - t89);
t13 = (t28 - t73) * t52;
t9 = (-t27 + t72) * t50;
t5 = -t57 + t93;
t4 = t58 - t91;
t1 = t4 * t50 + t5 * t52;
t2 = [0, 0, 0, 0, 0, qJDD(1), t74, t85, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t69 - 0.2e1 * t80, t42 + 0.2e1 * t47 - t85, pkin(1) * t18 + qJ(2) * (t42 + t63), t13, -t65, t14, t9, -t64, 0, -t16 * t50 + t71, qJ(2) * t29 - t52 * t16 - t98, -t3 + t70, -qJ(2) * t16 + t3 * t94, t13, t14, t65, 0, t64, t9, -t26 * t83 + t71 + (pkin(3) * t26 - t55 * pkin(5) + t63 - t95) * t50, t52 * (qJ(4) * t32 + t57) - t50 * (pkin(3) * t32 + t60 - t91) + t70, t98 - t62 * t29 + (t61 + t95) * t52, t94 * t1 + t62 * (qJD(1) * t99 + t42 - t61 - t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t55, -t18, 0, 0, 0, 0, 0, 0, t10, -t11, -t30, t3, 0, 0, 0, 0, 0, 0, t10, -t30, t11, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t31, t78, -t38, -t79, qJDD(3), t6, t7, 0, 0, t38, t78, -t31, qJDD(3), t79, -t38, t93 + (t36 + t54) * qJ(4) + (qJDD(3) + t34) * pkin(3) + t68, -t67 * qJDD(1), qJ(4) * t33 + (t37 - t54) * pkin(3) + t58, pkin(3) * t5 + qJ(4) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t78, -t37, -t5;];
tauJ_reg = t2;
