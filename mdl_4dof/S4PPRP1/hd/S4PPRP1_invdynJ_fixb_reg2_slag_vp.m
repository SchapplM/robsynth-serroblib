% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRP1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:16
% EndTime: 2018-11-14 13:39:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (80->41), mult. (125->48), div. (0->0), fcn. (98->4), ass. (0->26)
t18 = sin(qJ(3));
t28 = t18 * qJD(2);
t29 = qJDD(3) * pkin(3);
t31 = qJD(3) * t28 - t29;
t30 = sin(pkin(5));
t19 = cos(qJ(3));
t27 = t19 * qJD(2);
t26 = qJDD(3) * qJ(4);
t8 = qJD(3) * qJ(4) + t28;
t25 = t18 * (-qJD(3) * pkin(3) + qJD(4) - t27) + t19 * t8;
t17 = cos(pkin(5));
t24 = -g(1) * t30 + g(2) * t17;
t13 = t18 * qJDD(2);
t3 = -t17 * t19 - t30 * t18;
t4 = t17 * t18 - t30 * t19;
t23 = -g(1) * t3 - g(2) * t4 - t13;
t14 = t19 * qJDD(2);
t22 = g(1) * t4 - g(2) * t3 + t14;
t21 = -qJDD(4) + t22;
t20 = qJD(3) ^ 2;
t16 = qJDD(1) - g(3);
t6 = t19 * qJDD(3) - t20 * t18;
t5 = qJDD(3) * t18 + t20 * t19;
t2 = qJDD(4) - t14 + t31;
t1 = t26 + t13 + (qJD(4) + t27) * qJD(3);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t24, 0, 0, 0, 0, 0, 0, t6, -t5, 0 (t18 ^ 2 + t19 ^ 2) * qJDD(2) + t24, 0, 0, 0, 0, 0, 0, t6, 0, t5, t25 * qJD(3) + t1 * t18 - t2 * t19 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t22, t23, 0, 0, 0, 0, 0, qJDD(3), 0, 0, t21 + 0.2e1 * t29, 0, 0.2e1 * qJD(3) * qJD(4) - t23 + 0.2e1 * t26, t1 * qJ(4) + t8 * qJD(4) - t2 * pkin(3) - g(1) * (-t4 * pkin(3) - t3 * qJ(4)) - g(2) * (t3 * pkin(3) - t4 * qJ(4)) - t25 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3), 0, -t20, -t8 * qJD(3) - t21 + t31;];
tau_reg  = t7;
