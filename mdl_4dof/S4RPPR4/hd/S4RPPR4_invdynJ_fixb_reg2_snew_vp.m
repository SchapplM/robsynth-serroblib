% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (362->63), mult. (655->84), div. (0->0), fcn. (344->6), ass. (0->52)
t60 = -pkin(2) - pkin(5);
t35 = sin(qJ(4));
t29 = t35 ^ 2;
t40 = qJD(1) ^ 2;
t59 = t29 * t40;
t37 = cos(qJ(4));
t30 = t37 ^ 2;
t58 = t30 * t40;
t50 = t35 * t40 * t37;
t21 = qJDD(4) + t50;
t57 = t35 * t21;
t22 = qJDD(4) - t50;
t56 = t37 * t22;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t49 = t36 * g(1) - t38 * g(2);
t13 = qJDD(1) * pkin(1) + t49;
t44 = t38 * g(1) + t36 * g(2);
t14 = -t40 * pkin(1) - t44;
t33 = sin(pkin(6));
t34 = cos(pkin(6));
t55 = t33 * t13 + t34 * t14;
t54 = t29 + t30;
t53 = t35 * qJDD(1);
t52 = qJD(1) * qJD(4);
t51 = qJDD(1) * qJ(3);
t48 = pkin(1) * t33 + qJ(3);
t47 = t34 * t13 - t33 * t14;
t46 = pkin(1) * t34 - t60;
t45 = -pkin(1) * (t33 * qJDD(1) + t34 * t40) - t55;
t31 = -g(3) + qJDD(2);
t32 = qJDD(1) * pkin(2);
t7 = -t40 * qJ(3) + qJDD(3) - t32 - t47;
t41 = -qJDD(1) * pkin(5) + t7;
t2 = t35 * t31 - t37 * t41;
t3 = t37 * t31 + t35 * t41;
t1 = -t37 * t2 + t35 * t3;
t25 = t37 * qJDD(1);
t16 = -0.2e1 * t35 * t52 + t25;
t43 = -pkin(1) * (-t34 * qJDD(1) + t33 * t40) + t47;
t26 = 2 * qJD(3) * qJD(1);
t42 = t26 + t51 + t55;
t15 = 0.2e1 * t37 * t52 + t53;
t39 = qJD(4) ^ 2;
t24 = -t39 - t58;
t23 = -t39 - t59;
t19 = t54 * qJDD(1);
t9 = t37 * t24 - t57;
t8 = t35 * t23 + t56;
t6 = -t40 * pkin(2) + t42;
t5 = t60 * t40 + t42;
t4 = [0, 0, 0, 0, 0, qJDD(1), t49, t44, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t43, t45, 0, pkin(1) * (t33 * t55 + t34 * t47), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t32 - t43, t26 - t45 + 0.2e1 * t51, pkin(1) * (t33 * t6 - t34 * t7) - pkin(2) * t7 + qJ(3) * t6, t16 * t37, -t37 * t15 - t35 * t16, t56 - t35 * (t39 - t58), t15 * t35, t37 * (-t39 + t59) - t57, 0, t48 * t15 + t35 * t5 - t46 * t8, t48 * t16 + t37 * t5 - t46 * t9, -t48 * t54 * t40 + t46 * t19 - t1, -t46 * t1 + t48 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, -t35 * t22 + t37 * t23, -t37 * t21 - t35 * t24, 0, t35 * t2 + t37 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t40, t7, 0, 0, 0, 0, 0, 0, t8, t9, -t19, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, (-t29 + t30) * t40, t25, -t50, -t53, qJDD(4), -t2, -t3, 0, 0;];
tauJ_reg = t4;
