% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPPR5
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
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:22
% EndTime: 2018-11-14 14:05:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (94->33), mult. (173->49), div. (0->0), fcn. (154->4), ass. (0->23)
t15 = sin(pkin(5));
t25 = t15 ^ 2 * qJDD(1) - g(1);
t24 = qJD(1) * t15;
t23 = qJDD(1) * t15;
t16 = cos(pkin(5));
t22 = g(1) * t16 - g(2) * t15;
t10 = -t16 * qJDD(1) + qJDD(3);
t17 = sin(qJ(4));
t18 = cos(qJ(4));
t21 = t18 * t10 - t17 * t23;
t8 = t15 * t18 - t16 * t17;
t7 = -t15 * t17 - t16 * t18;
t11 = -t16 * qJD(1) + qJD(3);
t4 = t17 * t11 + t18 * t24;
t3 = t18 * t11 - t17 * t24;
t20 = t17 * t10 + t18 * t23;
t19 = qJD(4) ^ 2;
t14 = qJDD(2) + g(3);
t6 = t7 * qJD(4);
t5 = t8 * qJD(4);
t2 = t3 * qJD(4) + t20;
t1 = -t4 * qJD(4) + t21;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 * qJDD(1) + t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t16 + t25, 0, 0, 0, 0, 0, 0, -t5 * qJD(4) + t7 * qJDD(4), -t6 * qJD(4) - t8 * qJDD(4), 0, t1 * t7 + t2 * t8 - t3 * t5 + t4 * t6 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 + t22, 0, 0, 0, 0, 0, 0, t18 * qJDD(4) - t19 * t17, -qJDD(4) * t17 - t19 * t18, 0, t1 * t18 + t2 * t17 + (-t17 * t3 + t18 * t4) * qJD(4) + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -g(1) * t7 - g(2) * t8 + t21, g(1) * t8 - g(2) * t7 - t20, 0, 0;];
tau_reg  = t9;
