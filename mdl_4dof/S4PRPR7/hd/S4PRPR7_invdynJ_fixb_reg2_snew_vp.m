% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:56
% DurationCPUTime: 0.18s
% Computational Cost: add. (267->59), mult. (480->79), div. (0->0), fcn. (282->6), ass. (0->46)
t53 = pkin(2) + pkin(5);
t34 = sin(qJ(4));
t28 = t34 ^ 2;
t39 = qJD(2) ^ 2;
t52 = t28 * t39;
t36 = cos(qJ(4));
t29 = t36 ^ 2;
t51 = t29 * t39;
t44 = t36 * t39 * t34;
t50 = t34 * (qJDD(4) + t44);
t49 = t36 * (qJDD(4) - t44);
t32 = sin(pkin(6));
t33 = cos(pkin(6));
t20 = -t33 * g(1) - t32 * g(2);
t30 = -g(3) + qJDD(1);
t35 = sin(qJ(2));
t37 = cos(qJ(2));
t11 = t37 * t20 + t35 * t30;
t48 = t28 + t29;
t47 = t34 * qJDD(2);
t46 = qJD(2) * qJD(4);
t45 = qJDD(2) * qJ(3);
t43 = (2 * qJD(3) * qJD(2)) + t11;
t10 = -t35 * t20 + t37 * t30;
t19 = -t32 * g(1) + t33 * g(2);
t31 = qJDD(2) * pkin(2);
t42 = qJDD(3) - t10;
t7 = -t39 * qJ(3) - t31 + t42;
t40 = -qJDD(2) * pkin(5) + t7;
t2 = t34 * t19 - t36 * t40;
t3 = t36 * t19 + t34 * t40;
t1 = -t36 * t2 + t34 * t3;
t25 = t36 * qJDD(2);
t14 = -0.2e1 * t34 * t46 + t25;
t41 = t43 + t45;
t13 = 0.2e1 * t36 * t46 + t47;
t38 = qJD(4) ^ 2;
t18 = t48 * t39;
t17 = -t37 * qJDD(2) + t35 * t39;
t16 = t35 * qJDD(2) + t37 * t39;
t15 = t48 * qJDD(2);
t9 = -t50 + t36 * (-t38 - t51);
t8 = t34 * (-t38 - t52) + t49;
t6 = -t39 * pkin(2) + t41;
t5 = -t53 * t39 + t41;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, -t17, -t16, 0, t37 * t10 + t35 * t11, 0, 0, 0, 0, 0, 0, 0, t17, t16, t35 * t6 - t37 * t7, 0, 0, 0, 0, 0, 0, t35 * t13 - t37 * t8, t35 * t14 - t37 * t9, t37 * t15 - t35 * t18, -t37 * t1 + t35 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t10, -t11, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t31 + t42, t43 + 0.2e1 * t45, -pkin(2) * t7 + qJ(3) * t6, t14 * t36, -t36 * t13 - t34 * t14, t49 - t34 * (t38 - t51), t13 * t34, t36 * (-t38 + t52) - t50, 0, qJ(3) * t13 + t34 * t5 - t53 * t8, qJ(3) * t14 + t36 * t5 - t53 * t9, -qJ(3) * t18 + t53 * t15 - t1, qJ(3) * t5 - t53 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t39, t7, 0, 0, 0, 0, 0, 0, t8, t9, -t15, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, (-t28 + t29) * t39, t25, -t44, -t47, qJDD(4), -t2, -t3, 0, 0;];
tauJ_reg = t4;
