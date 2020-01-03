% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRP6
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:15
% DurationCPUTime: 0.36s
% Computational Cost: add. (445->103), mult. (938->108), div. (0->0), fcn. (419->4), ass. (0->78)
t78 = (qJD(1) * qJD(4));
t94 = 2 * t78;
t92 = -pkin(5) - pkin(1);
t77 = (qJD(2) * qJD(1));
t43 = 2 * t77;
t50 = sin(qJ(3));
t47 = t50 ^ 2;
t55 = qJD(1) ^ 2;
t52 = cos(qJ(3));
t79 = qJD(1) * qJD(3);
t70 = t52 * t79;
t80 = t50 * qJDD(1);
t29 = -t70 - t80;
t82 = qJD(1) * t52;
t35 = qJD(3) * pkin(3) - qJ(4) * t82;
t46 = qJDD(1) * qJ(2);
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t66 = t53 * g(1) + t51 * g(2);
t62 = -t46 + t66;
t58 = -t29 * pkin(3) + t35 * t82 + qJDD(4) - t62;
t5 = (-qJ(4) * t47 + t92) * t55 + t43 + t58;
t73 = t51 * g(1) - t53 * g(2);
t65 = qJDD(2) - t73;
t60 = -t55 * qJ(2) + t65;
t20 = t92 * qJDD(1) + t60;
t17 = t50 * t20;
t93 = -t29 * qJ(4) + qJD(3) * t35 + t50 * t94 - t17;
t44 = t50 * g(3);
t91 = t52 * g(3);
t42 = t52 * qJDD(1);
t71 = t50 * t79;
t30 = t42 - t71;
t63 = -t30 - t71;
t7 = t52 * t20 + t44;
t75 = qJDD(3) * pkin(3) + t7;
t57 = t63 * qJ(4) - 0.2e1 * t52 * t78 + t75;
t85 = t50 * t55;
t76 = t52 * t85;
t2 = -pkin(3) * t76 + t57;
t90 = t52 * t2;
t89 = t55 * pkin(1);
t88 = t47 * t55;
t48 = t52 ^ 2;
t87 = t48 * t55;
t36 = qJDD(3) + t76;
t86 = t50 * t36;
t37 = qJDD(3) - t76;
t84 = t52 * t37;
t83 = t47 + t48;
t81 = qJDD(1) * pkin(1);
t74 = -2 * t77;
t54 = qJD(3) ^ 2;
t38 = -t54 - t88;
t12 = t50 * t38 + t84;
t28 = 0.2e1 * t70 + t80;
t69 = qJ(2) * t28 + t92 * t12;
t39 = -t54 - t87;
t13 = t52 * t39 - t86;
t31 = t42 - 0.2e1 * t71;
t68 = qJ(2) * t31 + t92 * t13;
t32 = t83 * qJDD(1);
t34 = t83 * t55;
t67 = -qJ(2) * t34 - t92 * t32;
t8 = -t17 + t91;
t4 = -t50 * t8 + t52 * t7;
t59 = -pkin(3) * t88 - t93;
t33 = (-t47 + t48) * t55;
t21 = -t60 + t81;
t19 = -t92 * t55 + t62 + t74;
t16 = t84 - t50 * (t54 - t87);
t15 = (t30 - t71) * t52;
t14 = t52 * (-t54 + t88) - t86;
t11 = (-t29 + t70) * t50;
t6 = -t52 * t28 - t50 * t31;
t3 = t59 - t91;
t1 = t50 * t3 + t90;
t9 = [0, 0, 0, 0, 0, qJDD(1), t73, t66, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t65 - 0.2e1 * t81, t43 + 0.2e1 * t46 - t66, pkin(1) * t21 + qJ(2) * (t43 - t62 - t89), t15, t6, t16, t11, t14, 0, -t50 * t19 + t69, -t52 * t19 + t68, -t4 + t67, -qJ(2) * t19 + t92 * t4, t15, t6, t16, t11, t14, 0, -t50 * (-pkin(3) * t28 + t55 * pkin(5) - t58 + t74 + t89) + (-t84 - t50 * (t38 + t88)) * qJ(4) + t69, -t50 * (-pkin(3) * t31 - qJ(4) * t36) + t68 + (-qJ(4) * t39 + t5) * t52, -t50 * (pkin(3) * t34 - qJ(4) * t80 + t59) + (t44 + (pkin(3) * t85 + t94) * t52 + (-t63 + t42) * qJ(4) - t75) * t52 + t67, -qJ(4) * t90 - t50 * (-pkin(3) * t5 + qJ(4) * t3) + qJ(2) * t5 + t92 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t55, -t21, 0, 0, 0, 0, 0, 0, t12, t13, -t32, t4, 0, 0, 0, 0, 0, 0, t12, t13, -t32, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t33, t42, -t76, -t80, qJDD(3), t7, t8, 0, 0, t76, t33, t42, -t76, -t80, qJDD(3), (t37 - t76) * pkin(3) + t57, t91 + (t39 + t88) * pkin(3) + t93, -pkin(3) * t42, pkin(3) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, -t34, t5;];
tauJ_reg = t9;
