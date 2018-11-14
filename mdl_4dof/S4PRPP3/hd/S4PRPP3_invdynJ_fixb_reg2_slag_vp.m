% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:26
% EndTime: 2018-11-14 14:01:26
% DurationCPUTime: 0.16s
% Computational Cost: add. (121->55), mult. (174->56), div. (0->0), fcn. (90->2), ass. (0->34)
t24 = -pkin(2) - pkin(3);
t30 = t24 * qJDD(2);
t22 = sin(qJ(2));
t13 = qJD(2) * qJ(3) + qJD(1) * t22;
t15 = t22 * qJDD(1);
t33 = qJDD(2) * qJ(3);
t23 = cos(qJ(2));
t35 = t23 * qJD(1);
t5 = t33 + t15 + (qJD(3) + t35) * qJD(2);
t41 = t5 * qJ(3) + t13 * qJD(3);
t40 = t13 * t23;
t39 = t23 * pkin(2) + t22 * qJ(3);
t38 = qJD(2) * t22;
t37 = qJDD(2) * pkin(2);
t36 = t13 * qJD(2);
t34 = t5 * t22 + t23 * t36 - g(2);
t16 = t23 * qJDD(1);
t32 = g(1) * t22 - g(2) * t23 + t16;
t14 = qJD(1) * t38;
t31 = qJDD(3) + t14 - t16;
t29 = qJD(3) - t35;
t28 = -qJDD(3) + t32;
t27 = g(1) * t23 + g(2) * t22 - t15;
t26 = t14 - t28 - t36;
t25 = qJD(2) ^ 2;
t18 = t23 * qJ(3);
t12 = -qJD(2) * pkin(2) + t29;
t11 = qJDD(2) * t23 - t22 * t25;
t10 = qJDD(2) * t22 + t23 * t25;
t8 = t24 * qJD(2) + t29;
t6 = t31 - t37;
t4 = t30 + t31;
t2 = 0.2e1 * qJD(2) * qJD(3) - t27 + 0.2e1 * t33;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(2), 0, 0, 0, 0, 0, 0, t11, -t10, 0, -g(2) + (t22 ^ 2 + t23 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t11, 0, t10, t12 * t38 - t23 * t6 + t34, 0, 0, 0, 0, 0, 0, t11, t10, 0, -t23 * t4 + t8 * t38 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t32, t27, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t28 + 0.2e1 * t37, 0, t2, -t6 * pkin(2) - g(1) * (-pkin(2) * t22 + t18) - g(2) * t39 + (-t12 * t22 - t40) * qJD(1) + t41, 0, 0, 0, 0, 0, qJDD(2), t28 - 0.2e1 * t30, t2, 0, t4 * t24 - g(1) * (t24 * t22 + t18) - g(2) * (pkin(3) * t23 + t39) + (-t22 * t8 - t40) * qJD(1) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t25, t26 - t37, 0, 0, 0, 0, 0, 0, -qJDD(2), -t25, 0, t30 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t1;
