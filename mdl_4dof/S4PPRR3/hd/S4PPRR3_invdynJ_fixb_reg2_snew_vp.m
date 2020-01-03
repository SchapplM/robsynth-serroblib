% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (184->53), mult. (352->73), div. (0->0), fcn. (234->6), ass. (0->39)
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t40 = qJD(3) ^ 2;
t22 = t35 * t40 * t37;
t18 = qJDD(4) + t22;
t44 = t35 * t18;
t19 = qJDD(4) - t22;
t43 = t37 * t19;
t42 = t35 * qJDD(3);
t41 = qJD(3) * qJD(4);
t30 = g(3) - qJDD(1);
t32 = sin(pkin(6));
t33 = cos(pkin(6));
t14 = -t32 * g(1) + t33 * g(2) + qJDD(2);
t17 = -t33 * g(1) - t32 * g(2);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t8 = t36 * t14 + t38 * t17;
t6 = -t40 * pkin(3) + qJDD(3) * pkin(5) + t8;
t3 = -t37 * t30 + t35 * t6;
t4 = t35 * t30 + t37 * t6;
t1 = t35 * t3 + t37 * t4;
t7 = t38 * t14 - t36 * t17;
t25 = t37 * qJDD(3);
t13 = -0.2e1 * t35 * t41 + t25;
t12 = 0.2e1 * t37 * t41 + t42;
t39 = qJD(4) ^ 2;
t29 = t37 ^ 2;
t28 = t35 ^ 2;
t27 = t29 * t40;
t26 = t28 * t40;
t21 = -t27 - t39;
t20 = -t26 - t39;
t16 = t26 + t27;
t15 = (t28 + t29) * qJDD(3);
t10 = -t35 * t20 - t43;
t9 = t37 * t21 - t44;
t5 = -qJDD(3) * pkin(3) - t40 * pkin(5) - t7;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, -t37 * t18 - t35 * t21, t35 * t19 - t37 * t20, 0, t37 * t3 - t35 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, t38 * qJDD(3) - t36 * t40, -t36 * qJDD(3) - t38 * t40, 0, t36 * t8 + t38 * t7, 0, 0, 0, 0, 0, 0, t38 * t13 + t36 * t9, t36 * t10 - t38 * t12, t36 * t15 + t38 * t16, t36 * t1 - t38 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t7, -t8, 0, 0, t12 * t35, t37 * t12 + t35 * t13, t44 + t37 * (-t26 + t39), t13 * t37, t35 * (t27 - t39) + t43, 0, pkin(3) * t13 + pkin(5) * t9 - t37 * t5, -pkin(3) * t12 + pkin(5) * t10 + t35 * t5, pkin(3) * t16 + pkin(5) * t15 + t1, -pkin(3) * t5 + pkin(5) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t26 - t27, t42, t22, t25, qJDD(4), -t3, -t4, 0, 0;];
tauJ_reg = t2;
