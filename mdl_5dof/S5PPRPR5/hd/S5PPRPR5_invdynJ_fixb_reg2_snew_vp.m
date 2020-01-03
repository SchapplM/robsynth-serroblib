% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (326->67), mult. (556->85), div. (0->0), fcn. (342->6), ass. (0->50)
t55 = pkin(3) + pkin(6);
t36 = sin(qJ(5));
t30 = t36 ^ 2;
t41 = qJD(3) ^ 2;
t54 = t30 * t41;
t38 = cos(qJ(5));
t31 = t38 ^ 2;
t53 = t31 * t41;
t46 = t38 * t41 * t36;
t23 = qJDD(5) + t46;
t52 = t36 * t23;
t24 = qJDD(5) - t46;
t51 = t38 * t24;
t34 = sin(pkin(7));
t35 = cos(pkin(7));
t17 = -t34 * g(1) + t35 * g(2) + qJDD(2);
t22 = -t35 * g(1) - t34 * g(2);
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t9 = t37 * t17 + t39 * t22;
t50 = t30 + t31;
t49 = t36 * qJDD(3);
t48 = qJD(3) * qJD(5);
t47 = qJDD(3) * qJ(4);
t45 = (2 * qJD(4) * qJD(3)) + t9;
t8 = t39 * t17 - t37 * t22;
t32 = g(3) - qJDD(1);
t33 = qJDD(3) * pkin(3);
t44 = qJDD(4) - t8;
t7 = -t41 * qJ(4) - t33 + t44;
t42 = -qJDD(3) * pkin(6) + t7;
t2 = t36 * t32 - t38 * t42;
t3 = t38 * t32 + t36 * t42;
t1 = -t38 * t2 + t36 * t3;
t27 = t38 * qJDD(3);
t16 = -0.2e1 * t36 * t48 + t27;
t43 = t45 + t47;
t15 = 0.2e1 * t38 * t48 + t49;
t40 = qJD(5) ^ 2;
t26 = -t40 - t53;
t25 = -t40 - t54;
t21 = t50 * t41;
t20 = -t39 * qJDD(3) + t37 * t41;
t19 = t37 * qJDD(3) + t39 * t41;
t18 = t50 * qJDD(3);
t11 = t38 * t26 - t52;
t10 = t36 * t25 + t51;
t6 = -t41 * pkin(3) + t43;
t5 = -t55 * t41 + t43;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, t36 * t24 - t38 * t25, t38 * t23 + t36 * t26, 0, -t36 * t2 - t38 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, -t20, -t19, 0, t37 * t9 + t39 * t8, 0, 0, 0, 0, 0, 0, 0, t20, t19, t37 * t6 - t39 * t7, 0, 0, 0, 0, 0, 0, -t39 * t10 + t37 * t15, -t39 * t11 + t37 * t16, t39 * t18 - t37 * t21, -t39 * t1 + t37 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t8, -t9, 0, 0, qJDD(3), 0, 0, 0, 0, 0, 0, -0.2e1 * t33 + t44, t45 + 0.2e1 * t47, -pkin(3) * t7 + qJ(4) * t6, t16 * t38, -t38 * t15 - t36 * t16, t51 - t36 * (t40 - t53), t15 * t36, t38 * (-t40 + t54) - t52, 0, qJ(4) * t15 - t55 * t10 + t36 * t5, qJ(4) * t16 - t55 * t11 + t38 * t5, -qJ(4) * t21 + t55 * t18 - t1, qJ(4) * t5 - t55 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t41, t7, 0, 0, 0, 0, 0, 0, t10, t11, -t18, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, (-t30 + t31) * t41, t27, -t46, -t49, qJDD(5), -t2, -t3, 0, 0;];
tauJ_reg = t4;
