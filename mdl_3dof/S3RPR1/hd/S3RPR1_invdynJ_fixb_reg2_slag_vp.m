% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% tau_reg [3x(3*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3RPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:30
% EndTime: 2018-11-14 10:14:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (198->59), mult. (299->70), div. (0->0), fcn. (150->4), ass. (0->37)
t42 = qJD(1) - qJD(3);
t32 = -pkin(1) - pkin(2);
t15 = t32 * qJD(1) + qJD(2);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t45 = qJ(2) * qJD(1);
t5 = t30 * t15 - t28 * t45;
t56 = t42 * t5;
t55 = -qJD(3) * t45 + t32 * qJDD(1) + qJDD(2);
t43 = qJD(1) * qJD(2);
t44 = qJ(2) * qJDD(1);
t54 = qJD(3) * t15 + t43 + t44;
t2 = -t54 * t28 + t55 * t30;
t6 = t28 * t15 + t30 * t45;
t53 = -t6 * t42 + t2;
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t8 = -t29 * t28 - t31 * t30;
t9 = t31 * t28 - t29 * t30;
t52 = g(1) * t9 - g(2) * t8;
t49 = t31 * pkin(1) + t29 * qJ(2);
t48 = g(1) * t29 - g(2) * t31;
t47 = pkin(1) * qJDD(1);
t41 = 0.2e1 * t43;
t39 = t42 ^ 2;
t38 = qJDD(2) - t47;
t37 = g(1) * t31 + g(2) * t29;
t13 = t30 * qJ(2) + t28 * t32;
t12 = -t28 * qJ(2) + t30 * t32;
t1 = t55 * t28 + t54 * t30;
t34 = -g(1) * t8 - g(2) * t9 - t1;
t33 = qJD(1) ^ 2;
t26 = qJDD(1) - qJDD(3);
t22 = t31 * qJ(2);
t4 = -t28 * qJD(2) - t13 * qJD(3);
t3 = t30 * qJD(2) + t12 * qJD(3);
t7 = [0, 0, 0, 0, 0, qJDD(1), t48, t37, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t47 + t48, 0, -t37 + t41 + 0.2e1 * t44, -t38 * pkin(1) - g(1) * (-t29 * pkin(1) + t22) - g(2) * t49 + (t41 + t44) * qJ(2), 0, 0, 0, 0, 0, t26, -t12 * t26 - t4 * t42 - t2 - t52, t13 * t26 + t3 * t42 - t34, 0, t1 * t13 + t6 * t3 + t2 * t12 + t5 * t4 - g(1) * (t32 * t29 + t22) - g(2) * (t31 * pkin(2) + t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t33, -t33 * qJ(2) + t38 - t48, 0, 0, 0, 0, 0, 0, -t30 * t26 - t28 * t39, t28 * t26 - t30 * t39, 0, t53 * t30 + (t1 + t56) * t28 - t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t52 + t53, t34 - t56, 0, 0;];
tau_reg  = t7;
