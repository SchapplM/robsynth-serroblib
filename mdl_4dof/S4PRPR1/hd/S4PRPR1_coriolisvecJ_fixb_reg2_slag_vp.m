% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:23
% EndTime: 2018-11-14 13:42:23
% DurationCPUTime: 0.11s
% Computational Cost: add. (109->26), mult. (192->39), div. (0->0), fcn. (74->2), ass. (0->25)
t12 = cos(qJ(4));
t11 = sin(qJ(4));
t24 = qJD(2) * qJ(3);
t19 = t11 * t24;
t23 = qJD(2) * qJD(3);
t26 = qJD(4) * t12;
t13 = -pkin(2) - pkin(3);
t8 = t13 * qJD(2) + qJD(3);
t1 = -qJD(4) * t19 + t12 * t23 + t8 * t26;
t22 = qJD(2) - qJD(4);
t5 = t12 * t8 - t19;
t29 = t22 * t5 + t1;
t25 = t11 * qJD(3);
t27 = t11 * t8;
t2 = -(qJ(3) * t26 + t25) * qJD(2) - qJD(4) * t27;
t6 = t12 * t24 + t27;
t28 = -t6 * t22 + t2;
t20 = 0.2e1 * t23;
t18 = t22 ^ 2;
t17 = t12 * qJ(3) + t11 * t13;
t16 = -t11 * qJ(3) + t12 * t13;
t14 = qJD(2) ^ 2;
t4 = -t17 * qJD(4) - t25;
t3 = t12 * qJD(3) + t16 * qJD(4);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, qJ(3) * t20, 0, 0, 0, 0, 0, 0, -t22 * t4 - t2, t22 * t3 + t1, 0, t1 * t17 + t2 * t16 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t14 * qJ(3), 0, 0, 0, 0, 0, 0, -t11 * t18, -t12 * t18, 0, t29 * t11 + t28 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, 0;];
tauc_reg  = t7;
