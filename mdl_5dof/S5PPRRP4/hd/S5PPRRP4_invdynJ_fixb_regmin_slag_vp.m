% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:43
% EndTime: 2019-12-31 17:34:44
% DurationCPUTime: 0.39s
% Computational Cost: add. (365->109), mult. (763->149), div. (0->0), fcn. (484->6), ass. (0->66)
t72 = qJ(5) + pkin(6);
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t38 = qJD(3) ^ 2;
t76 = qJDD(3) * t34 + t38 * t36;
t35 = cos(qJ(4));
t26 = pkin(4) * t35 + pkin(3);
t53 = qJD(2) * qJD(3);
t24 = t34 * t53;
t56 = t36 * qJDD(2);
t48 = t24 - t56;
t33 = sin(qJ(4));
t52 = qJD(3) * qJD(4);
t51 = t33 * t52;
t3 = pkin(4) * t51 - qJDD(3) * t26 + qJDD(5) + t48;
t49 = t34 * qJD(2) + qJD(3) * t72;
t5 = -t35 * qJD(1) - t33 * t49;
t66 = qJD(4) * pkin(4);
t4 = t5 + t66;
t27 = t33 * qJD(1);
t6 = t35 * t49 - t27;
t45 = t33 * t4 - t35 * t6;
t75 = -qJD(3) * t45 - t3;
t74 = t4 - t5;
t29 = t33 ^ 2;
t30 = t35 ^ 2;
t71 = t29 - t30;
t70 = t29 + t30;
t37 = qJD(4) ^ 2;
t69 = t37 + t38;
t68 = qJD(3) * pkin(3);
t62 = t36 * qJD(2);
t9 = -qJD(3) * t26 + qJD(5) - t62;
t67 = qJD(3) * t9;
t65 = cos(pkin(7));
t64 = sin(pkin(7));
t31 = qJDD(1) - g(3);
t60 = qJDD(4) * t33;
t59 = t33 * qJDD(3);
t58 = t34 * qJDD(2);
t57 = t35 * qJDD(3);
t55 = t36 * qJDD(3);
t54 = qJD(1) * qJD(4);
t50 = qJD(4) * t72;
t12 = -t34 * t64 - t36 * t65;
t13 = t34 * t65 - t36 * t64;
t47 = g(1) * t13 - g(2) * t12;
t46 = g(1) * t12 + g(2) * t13;
t44 = -g(1) * t64 + g(2) * t65;
t11 = qJDD(3) * pkin(6) + t36 * t53 + t58;
t43 = -qJ(5) * qJDD(3) - qJD(3) * qJD(5) - t11;
t42 = -qJD(4) * t49 - qJDD(1);
t19 = -t62 - t68;
t41 = -qJD(3) * t19 - t11 - t46;
t40 = -pkin(6) * qJDD(4) + (t19 + t62 - t68) * qJD(4);
t39 = 0.2e1 * qJDD(3) * pkin(3) - pkin(6) * t37 + t24 + t47 - t48;
t25 = t33 * t54;
t17 = t72 * t35;
t16 = t72 * t33;
t15 = qJDD(4) * t35 - t33 * t37;
t14 = t35 * t37 + t60;
t8 = -t33 * qJD(5) - t35 * t50;
t7 = t35 * qJD(5) - t33 * t50;
t2 = t42 * t33 + (-t43 - t54) * t35;
t1 = qJDD(4) * pkin(4) + t33 * t43 + t35 * t42 + t25;
t10 = [t31, t31, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t14, 0, qJD(4) * t45 - t1 * t35 - t2 * t33 - g(3); 0, qJDD(2) + t44, 0, -t34 * t38 + t55, -t76, 0, 0, 0, 0, 0, (-0.2e1 * t51 + t57) * t36 + (-t35 * t69 - t60) * t34, (-qJDD(4) * t34 - 0.2e1 * t36 * t52) * t35 + (t34 * t69 - t55) * t33, t76 * t70, t75 * t36 + (t67 - t1 * t33 + t2 * t35 + (-t33 * t6 - t35 * t4) * qJD(4)) * t34 + t44; 0, 0, qJDD(3), t47 + t56, -t46 - t58, qJDD(3) * t29 + 0.2e1 * t35 * t51, 0.2e1 * t33 * t57 - 0.2e1 * t52 * t71, t14, t15, 0, t33 * t40 + t35 * t39, -t33 * t39 + t35 * t40, (-qJD(4) * t4 + qJDD(3) * t17 + t2) * t35 + (-qJD(4) * t6 + qJDD(3) * t16 - t1) * t33 + (-t33 * t8 + t35 * t7 + (t16 * t35 - t17 * t33) * qJD(4) - t70 * t62) * qJD(3) + t46, t2 * t17 + t6 * t7 - t1 * t16 + t4 * t8 - t3 * t26 + t9 * t33 * t66 - g(1) * (-t12 * t72 - t13 * t26) - g(2) * (t12 * t26 - t13 * t72) + (-t34 * t9 + t36 * t45) * qJD(2); 0, 0, 0, 0, 0, -t33 * t38 * t35, t71 * t38, t59, t57, qJDD(4), -t27 * qJD(4) - t31 * t35 + t33 * t41 + t25, t31 * t33 + t35 * t41, -pkin(4) * t59 + (-t66 + t74) * t35 * qJD(3), t74 * t6 + (g(3) * t35 + t1 + (-t46 - t67) * t33) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 * t38, -t47 - t75;];
tau_reg = t10;
