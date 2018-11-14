% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRR1
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
% tau_reg [4x8]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:20
% EndTime: 2018-11-14 13:40:20
% DurationCPUTime: 0.15s
% Computational Cost: add. (112->47), mult. (170->64), div. (0->0), fcn. (144->8), ass. (0->37)
t13 = qJD(3) + qJD(4);
t39 = t13 ^ 2;
t15 = qJ(3) + qJ(4);
t38 = cos(t15);
t12 = qJDD(3) + qJDD(4);
t37 = pkin(3) * t12;
t18 = sin(qJ(4));
t21 = cos(qJ(3));
t36 = t18 * t21;
t19 = sin(qJ(3));
t35 = qJD(2) * t19;
t34 = qJD(4) * t19;
t33 = t21 * qJD(2);
t32 = -qJD(4) + t13;
t31 = t19 * qJDD(2);
t30 = qJD(2) * qJD(3);
t10 = t21 * qJDD(2);
t4 = qJDD(3) * pkin(3) - t19 * t30 + t10;
t29 = -t13 * t35 - t4;
t20 = cos(qJ(4));
t28 = t20 * t19 + t36;
t27 = -t18 * t19 + t20 * t21;
t11 = sin(t15);
t16 = sin(pkin(6));
t17 = cos(pkin(6));
t2 = -t16 * t11 - t17 * t38;
t3 = t17 * t11 - t16 * t38;
t26 = g(1) * t3 - g(2) * t2 + t20 * t4;
t25 = t18 * qJD(2) * t34 - g(1) * t2 - g(2) * t3;
t8 = qJD(3) * pkin(3) + t33;
t24 = -t31 + (-pkin(3) * t13 - t8) * qJD(4);
t23 = -t21 * t30 + t32 * t8 - t31;
t22 = qJD(3) ^ 2;
t14 = qJDD(1) - g(3);
t6 = -t16 * t21 + t17 * t19;
t5 = -t16 * t19 - t17 * t21;
t1 = [t14, t14, 0, 0, 0, 0, 0, 0; 0, -g(1) * t16 + g(2) * t17 + qJDD(2), 0, t21 * qJDD(3) - t22 * t19, -qJDD(3) * t19 - t22 * t21, 0, t27 * t12 - t39 * t28, -t28 * t12 - t39 * t27; 0, 0, qJDD(3), g(1) * t6 - g(2) * t5 + t10, -g(1) * t5 - g(2) * t6 - t31, t12, t20 * t37 + t24 * t18 + (-qJD(3) * t36 + t28 * t13 - t20 * t34) * qJD(2) + t26 (t29 - t37) * t18 + ((-qJD(3) + t13) * t33 + t24) * t20 + t25; 0, 0, 0, 0, 0, t12, t32 * t20 * t35 + t23 * t18 + t26, t29 * t18 + t23 * t20 + t25;];
tau_reg  = t1;
