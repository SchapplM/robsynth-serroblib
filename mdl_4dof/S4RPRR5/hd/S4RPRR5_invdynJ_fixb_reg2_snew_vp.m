% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:41
% DurationCPUTime: 0.29s
% Computational Cost: add. (875->91), mult. (1269->102), div. (0->0), fcn. (534->6), ass. (0->63)
t77 = pkin(1) + pkin(2);
t47 = (-qJD(1) + qJD(3));
t45 = t47 ^ 2;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t34 = t54 * t45 * t51;
t76 = t51 * (qJDD(4) + t34);
t46 = qJDD(1) - qJDD(3);
t75 = t51 * t46;
t74 = t54 * (qJDD(4) - t34);
t58 = qJD(1) ^ 2;
t48 = qJDD(1) * qJ(2);
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t68 = t56 * g(1) + t53 * g(2);
t65 = (2 * qJD(2) * qJD(1)) - t68;
t64 = t48 + t65;
t24 = -t77 * t58 + t64;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t71 = t53 * g(1) - t56 * g(2);
t66 = qJDD(2) - t71;
t61 = -t58 * qJ(2) + t66;
t59 = -t77 * qJDD(1) + t61;
t12 = t55 * t24 + t52 * t59;
t73 = qJDD(1) * pkin(1);
t72 = 2 * qJD(4) * t47;
t8 = -t45 * pkin(3) - t46 * pkin(6) + t12;
t5 = -t54 * g(3) + t51 * t8;
t6 = t51 * g(3) + t54 * t8;
t2 = t51 * t5 + t54 * t6;
t70 = t52 * t24 - t55 * t59;
t7 = t46 * pkin(3) - t45 * pkin(6) + t70;
t69 = -pkin(3) * t7 + pkin(6) * t2;
t36 = t54 * t46;
t67 = t51 * t72 + t36;
t29 = -t55 * t45 + t52 * t46;
t30 = t52 * t45 + t55 * t46;
t26 = t54 * t72 - t75;
t50 = t54 ^ 2;
t39 = t50 * t45;
t57 = qJD(4) ^ 2;
t19 = t54 * (-t39 - t57) - t76;
t63 = -pkin(3) * t67 + pkin(6) * t19 - t54 * t7;
t49 = t51 ^ 2;
t38 = t49 * t45;
t20 = -t74 - t51 * (-t38 - t57);
t62 = pkin(3) * t26 - pkin(6) * t20 - t51 * t7;
t28 = (-t49 - t50) * t46;
t31 = t38 + t39;
t60 = pkin(3) * t31 + pkin(6) * t28 + t2;
t25 = -t61 + t73;
t18 = t76 + t54 * (-t38 + t57);
t17 = t51 * (t39 - t57) + t74;
t16 = t26 * t51;
t15 = t67 * t54;
t14 = t52 * t28 + t55 * t31;
t13 = t54 * t26 - t51 * t67;
t10 = t52 * t20 - t55 * t26;
t9 = t52 * t19 - t55 * t67;
t3 = t52 * t12 - t55 * t70;
t1 = t52 * t2 - t55 * t7;
t4 = [0, 0, 0, 0, 0, qJDD(1), t71, t68, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t66 + 0.2e1 * t73, 0, 0.2e1 * t48 + t65, qJ(2) * (-t58 * pkin(1) + t64) + pkin(1) * t25, 0, 0, 0, 0, 0, t46, qJ(2) * t29 + t77 * t30 + t70, qJ(2) * t30 - t77 * t29 + t12, 0, qJ(2) * (t55 * t12 + t52 * t70) - t77 * t3, -t16, -t13, -t18, t15, -t17, 0, qJ(2) * (t55 * t19 + t52 * t67) - t77 * t9 - t63, qJ(2) * (t55 * t20 + t52 * t26) - t77 * t10 + t62, qJ(2) * (t55 * t28 - t52 * t31) - t77 * t14 - t60, qJ(2) * (t55 * t2 + t52 * t7) - t77 * t1 - t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t58, -t25, 0, 0, 0, 0, 0, 0, -t30, t29, 0, t3, 0, 0, 0, 0, 0, 0, t9, t10, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t70, -t12, 0, 0, t16, t13, t18, -t15, t17, 0, t63, -t62, t60, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t38 - t39, -t75, t34, -t36, qJDD(4), -t5, -t6, 0, 0;];
tauJ_reg = t4;
