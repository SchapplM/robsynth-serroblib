% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRPR2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:23
% DurationCPUTime: 0.30s
% Computational Cost: add. (539->80), mult. (897->119), div. (0->0), fcn. (638->8), ass. (0->56)
t62 = pkin(3) + pkin(6);
t41 = sin(qJ(5));
t33 = t41 ^ 2;
t46 = qJD(3) ^ 2;
t61 = t33 * t46;
t43 = cos(qJ(5));
t34 = t43 ^ 2;
t60 = t34 * t46;
t53 = t43 * t46 * t41;
t25 = qJDD(5) + t53;
t59 = t41 * t25;
t26 = qJDD(5) - t53;
t58 = t43 * t26;
t38 = sin(pkin(7));
t40 = cos(pkin(7));
t24 = -t40 * g(1) - t38 * g(2);
t35 = -g(3) + qJDD(1);
t37 = sin(pkin(8));
t39 = cos(pkin(8));
t15 = -t37 * t24 + t39 * t35;
t16 = t39 * t24 + t37 * t35;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t9 = t42 * t15 + t44 * t16;
t57 = t33 + t34;
t56 = t41 * qJDD(3);
t55 = qJD(3) * qJD(5);
t54 = qJDD(3) * qJ(4);
t52 = (2 * qJD(4) * qJD(3)) + t9;
t8 = t44 * t15 - t42 * t16;
t19 = -t38 * g(1) + t40 * g(2) + qJDD(2);
t36 = qJDD(3) * pkin(3);
t51 = qJDD(4) - t8;
t7 = -t46 * qJ(4) - t36 + t51;
t47 = -qJDD(3) * pkin(6) + t7;
t2 = t41 * t19 - t43 * t47;
t3 = t43 * t19 + t41 * t47;
t1 = -t43 * t2 + t41 * t3;
t29 = t43 * qJDD(3);
t18 = -0.2e1 * t41 * t55 + t29;
t21 = t44 * qJDD(3) - t42 * t46;
t22 = t42 * qJDD(3) + t44 * t46;
t50 = t39 * t21 - t37 * t22;
t49 = t37 * t21 + t39 * t22;
t48 = t52 + t54;
t17 = 0.2e1 * t43 * t55 + t56;
t45 = qJD(5) ^ 2;
t28 = -t45 - t60;
t27 = -t45 - t61;
t23 = t57 * t46;
t20 = t57 * qJDD(3);
t11 = t43 * t28 - t59;
t10 = t41 * t27 + t58;
t6 = -t46 * pkin(3) + t48;
t5 = -t62 * t46 + t48;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t15 + t37 * t16, 0, 0, 0, 0, 0, 0, t50, -t49, 0, t37 * (-t42 * t8 + t44 * t9) + t39 * (t42 * t9 + t44 * t8), 0, 0, 0, 0, 0, 0, 0, -t50, t49, t37 * (t42 * t7 + t44 * t6) + t39 * (t42 * t6 - t44 * t7), 0, 0, 0, 0, 0, 0, t37 * (t42 * t10 + t44 * t17) + t39 * (-t44 * t10 + t42 * t17), t37 * (t42 * t11 + t44 * t18) + t39 * (-t44 * t11 + t42 * t18), t37 * (-t42 * t20 - t44 * t23) + t39 * (t44 * t20 - t42 * t23), t37 * (t42 * t1 + t44 * t5) + t39 * (-t44 * t1 + t42 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, -t41 * t26 + t43 * t27, -t43 * t25 - t41 * t28, 0, t41 * t2 + t43 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t8, -t9, 0, 0, qJDD(3), 0, 0, 0, 0, 0, 0, -0.2e1 * t36 + t51, t52 + 0.2e1 * t54, -pkin(3) * t7 + qJ(4) * t6, t18 * t43, -t43 * t17 - t41 * t18, t58 - t41 * (t45 - t60), t17 * t41, t43 * (-t45 + t61) - t59, 0, qJ(4) * t17 - t62 * t10 + t41 * t5, qJ(4) * t18 - t62 * t11 + t43 * t5, -qJ(4) * t23 + t62 * t20 - t1, qJ(4) * t5 - t62 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t46, t7, 0, 0, 0, 0, 0, 0, t10, t11, -t20, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, (-t33 + t34) * t46, t29, -t53, -t56, qJDD(5), -t2, -t3, 0, 0;];
tauJ_reg = t4;
