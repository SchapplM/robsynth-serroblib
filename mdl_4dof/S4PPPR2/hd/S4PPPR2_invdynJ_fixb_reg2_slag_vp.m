% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPPR2
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
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:56
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:55:37
% EndTime: 2018-11-14 13:55:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (94->33), mult. (173->49), div. (0->0), fcn. (154->4), ass. (0->23)
t16 = sin(pkin(5));
t26 = t16 ^ 2 * qJDD(1) - g(2);
t25 = qJD(1) * t16;
t24 = qJDD(1) * t16;
t17 = cos(pkin(5));
t23 = -g(1) * t16 + g(2) * t17;
t10 = -t17 * qJDD(1) + qJDD(3);
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t22 = t19 * t10 - t18 * t24;
t8 = -t16 * t19 + t17 * t18;
t7 = -t16 * t18 - t17 * t19;
t11 = -t17 * qJD(1) + qJD(3);
t4 = t18 * t11 + t19 * t25;
t3 = t19 * t11 - t18 * t25;
t21 = t18 * t10 + t19 * t24;
t20 = qJD(4) ^ 2;
t15 = qJDD(2) - g(3);
t6 = t7 * qJD(4);
t5 = t8 * qJD(4);
t2 = t3 * qJD(4) + t21;
t1 = -t4 * qJD(4) + t22;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 * qJDD(1) + t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t17 + t26, 0, 0, 0, 0, 0, 0, t5 * qJD(4) + t7 * qJDD(4), -t6 * qJD(4) + t8 * qJDD(4), 0, t1 * t7 - t2 * t8 + t3 * t5 + t4 * t6 - g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 + t23, 0, 0, 0, 0, 0, 0, t19 * qJDD(4) - t20 * t18, -qJDD(4) * t18 - t20 * t19, 0, t1 * t19 + t2 * t18 + (-t18 * t3 + t19 * t4) * qJD(4) + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), g(1) * t8 - g(2) * t7 + t22, -g(1) * t7 - g(2) * t8 - t21, 0, 0;];
tau_reg  = t9;
