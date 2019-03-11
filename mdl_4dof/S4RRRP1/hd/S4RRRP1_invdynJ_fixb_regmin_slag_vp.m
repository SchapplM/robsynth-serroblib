% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:04
% EndTime: 2019-03-08 18:36:05
% DurationCPUTime: 0.31s
% Computational Cost: add. (339->95), mult. (547->118), div. (0->0), fcn. (291->10), ass. (0->67)
t35 = qJD(1) + qJD(2);
t41 = cos(qJ(2));
t67 = pkin(1) * qJD(1);
t14 = t35 * pkin(2) + t41 * t67;
t37 = sin(qJ(3));
t38 = sin(qJ(2));
t63 = qJD(3) * t38;
t55 = qJD(1) * t63;
t52 = pkin(1) * t55;
t15 = t37 * t52;
t40 = cos(qJ(3));
t66 = qJD(2) * t41;
t56 = qJD(1) * t66;
t61 = qJDD(1) * t38;
t46 = (t56 + t61) * pkin(1);
t72 = t41 * pkin(1);
t28 = qJDD(1) * t72;
t34 = qJDD(1) + qJDD(2);
t60 = t38 * t67;
t9 = t34 * pkin(2) - qJD(2) * t60 + t28;
t2 = t37 * t9 + (qJD(3) * t14 + t46) * t40 - t15;
t36 = qJ(1) + qJ(2);
t33 = qJ(3) + t36;
t25 = sin(t33);
t26 = cos(t33);
t76 = g(1) * t26 + g(2) * t25;
t18 = g(1) * t25;
t73 = g(2) * t26;
t75 = t18 - t73;
t29 = qJDD(3) + t34;
t74 = pkin(2) * t29;
t31 = sin(t36);
t22 = g(1) * t31;
t71 = t37 * t38;
t70 = t38 * t40;
t32 = cos(t36);
t69 = pkin(2) * t32 + pkin(3) * t26;
t68 = g(1) * t32 + g(2) * t31;
t64 = qJD(3) * t37;
t62 = qJD(3) * t40;
t59 = t38 * t62;
t8 = t40 * t9;
t50 = -t14 * t64 + t8;
t43 = (-t37 * t61 + (-t37 * t66 - t59) * qJD(1)) * pkin(1) + t50;
t1 = t29 * pkin(3) + t43;
t58 = t1 + t18;
t54 = qJD(1) * (-qJD(2) + t35);
t53 = qJD(2) * (-qJD(1) - t35);
t51 = -g(2) * t32 + t22 + t28;
t49 = -t37 * t41 - t70;
t48 = t40 * t41 - t71;
t6 = t40 * t14 - t37 * t60;
t30 = qJD(3) + t35;
t45 = (-pkin(2) * t30 - t14) * qJD(3) - t46;
t44 = -t2 + t76;
t42 = cos(qJ(1));
t39 = sin(qJ(1));
t27 = pkin(2) + t72;
t16 = t40 * t27;
t12 = pkin(1) * t70 + t37 * t27;
t11 = t48 * t67;
t10 = t49 * t67;
t7 = t37 * t14 + t40 * t60;
t5 = t30 * pkin(3) + t6;
t4 = -t27 * t64 + (t49 * qJD(2) - t59) * pkin(1);
t3 = t27 * t62 + (t48 * qJD(2) - t37 * t63) * pkin(1);
t13 = [qJDD(1), g(1) * t39 - g(2) * t42, g(1) * t42 + g(2) * t39, t34 (t34 * t41 + t38 * t53) * pkin(1) + t51 ((-qJDD(1) - t34) * t38 + t41 * t53) * pkin(1) + t68, t29, t16 * t29 + t4 * t30 + (-t40 * t55 + (-t56 + (-qJDD(1) - t29) * t38) * t37) * pkin(1) + t50 + t75, -t12 * t29 - t3 * t30 + t44, t2 * t12 + t7 * t3 + t1 * (-pkin(1) * t71 + pkin(3) + t16) + t5 * t4 - g(1) * (-t39 * pkin(1) - pkin(2) * t31 - pkin(3) * t25) - g(2) * (t42 * pkin(1) + t69); 0, 0, 0, t34, t38 * pkin(1) * t54 + t51 (t41 * t54 - t61) * pkin(1) + t68, t29, -t10 * t30 + t8 + (-t52 + t74) * t40 + t45 * t37 + t75, t11 * t30 + t15 + (-t9 - t74) * t37 + t45 * t40 + t76, -t7 * t11 - t5 * t10 - g(2) * t69 + t58 * pkin(3) + (t22 + t1 * t40 + t2 * t37 + (-t37 * t5 + t40 * t7) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, t29, t7 * t30 + t43 + t75, t6 * t30 + t44 (t5 - t6) * t7 + (t58 - t73) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) - g(3);];
tau_reg  = t13;
