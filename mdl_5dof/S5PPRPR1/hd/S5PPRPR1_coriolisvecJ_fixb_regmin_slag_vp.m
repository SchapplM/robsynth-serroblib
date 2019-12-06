% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:20
% EndTime: 2019-12-05 15:01:22
% DurationCPUTime: 0.32s
% Computational Cost: add. (273->60), mult. (799->100), div. (0->0), fcn. (626->8), ass. (0->45)
t40 = sin(pkin(8));
t42 = cos(pkin(8));
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t69 = -t44 * t40 + t46 * t42;
t70 = t69 * qJD(1);
t26 = t46 * t40 + t44 * t42;
t22 = t26 * qJD(3);
t39 = sin(pkin(9));
t41 = cos(pkin(9));
t59 = t39 ^ 2 + t41 ^ 2;
t12 = qJD(1) * t22;
t18 = t26 * qJD(1);
t50 = t18 * qJD(3) - t12;
t45 = cos(qJ(5));
t63 = t45 * t41;
t43 = sin(qJ(5));
t66 = t43 * t39;
t23 = -t63 + t66;
t19 = t23 * qJD(5);
t25 = t45 * t39 + t43 * t41;
t17 = t25 * qJD(3);
t67 = t70 - qJD(4);
t60 = pkin(6) + qJ(4);
t57 = t19 * qJD(5);
t21 = t69 * qJD(3);
t56 = t21 * qJD(3);
t55 = t22 * qJD(3);
t53 = qJD(3) * t66;
t52 = qJD(3) * t63;
t34 = -t41 * pkin(4) - pkin(3);
t7 = (qJD(4) + t70) * qJD(3);
t51 = t59 * t7;
t49 = t59 * (qJD(3) * qJ(4) + t18);
t20 = t25 * qJD(5);
t29 = qJD(5) * t52;
t28 = t60 * t41;
t27 = t60 * t39;
t14 = -t52 + t53;
t13 = t20 * qJD(5);
t11 = qJD(3) * t20;
t10 = -qJD(5) * t53 + t29;
t8 = -qJD(3) * pkin(3) - t67;
t5 = t34 * qJD(3) - t67;
t1 = [0, 0, 0, -t55, -t56, -t41 * t55, t39 * t55, t59 * t56, -t12 * t69 + t49 * t21 + t8 * t22 + t26 * t51, 0, 0, 0, 0, 0, -t69 * t11 + t22 * t14 + (t26 * t19 - t25 * t21) * qJD(5), -t69 * t10 + t22 * t17 + (t26 * t20 + t23 * t21) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t57; 0, 0, 0, t50, 0, t50 * t41, -t50 * t39, -t67 * qJD(3) * t59 + t51, -t12 * pkin(3) + qJ(4) * t51 - t8 * t18 - t67 * t49, t10 * t25 - t17 * t19, -t10 * t23 - t25 * t11 + t19 * t14 - t17 * t20, -t57, -t13, 0, t34 * t11 + t12 * t23 - t18 * t14 + t5 * t20 + ((t27 * t43 - t28 * t45) * qJD(5) + t67 * t25) * qJD(5), t34 * t10 + t12 * t25 - t18 * t17 - t5 * t19 + ((t27 * t45 + t28 * t43) * qJD(5) - t67 * t23) * qJD(5); 0, 0, 0, 0, 0, 0, 0, -t59 * qJD(3) ^ 2, -t49 * qJD(3) + t12, 0, 0, 0, 0, 0, 0.2e1 * t17 * qJD(5), t29 + (-t14 - t53) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t14, -t14 ^ 2 + t17 ^ 2, t29 + (t14 - t53) * qJD(5), 0, 0, -t5 * t17 - t25 * t7, t5 * t14 + t23 * t7;];
tauc_reg = t1;
