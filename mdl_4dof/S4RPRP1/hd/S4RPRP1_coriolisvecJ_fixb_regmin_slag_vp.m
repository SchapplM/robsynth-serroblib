% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:39
% EndTime: 2018-11-14 13:48:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (119->30), mult. (275->41), div. (0->0), fcn. (132->4), ass. (0->28)
t34 = pkin(1) * sin(pkin(6));
t21 = cos(qJ(3));
t28 = qJD(3) * t34;
t25 = qJD(1) * t28;
t15 = cos(pkin(6)) * pkin(1) + pkin(2);
t14 = t15 * qJD(1);
t20 = sin(qJ(3));
t33 = t20 * t14;
t4 = qJD(3) * t33 + t21 * t25;
t30 = qJD(3) * t21;
t32 = -t14 * t30 + t20 * t25;
t29 = qJD(1) * t34;
t5 = t21 * t14 - t20 * t29;
t31 = qJD(4) - t5;
t17 = qJD(1) + qJD(3);
t16 = t17 * qJD(4);
t1 = t16 - t32;
t6 = t21 * t29 + t33;
t27 = t6 * t17 - t4;
t22 = t20 * t15 + t21 * t34;
t8 = t22 * qJD(3);
t26 = -t8 * t17 - t4;
t24 = t5 * t17 + t32;
t23 = t15 * t30 - t20 * t28;
t7 = qJD(4) + t23;
t3 = t17 * qJ(4) + t6;
t2 = -t17 * pkin(3) + t31;
t9 = [0, 0, 0, 0, 0, t26, -t23 * t17 + t32, t26, t7 * t17 + t1, t1 * (qJ(4) + t22) + t3 * t7 + t4 * (-t21 * t15 + t20 * t34 - pkin(3)) + t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t27, t24, t27, 0.2e1 * t16 - t24, -t4 * pkin(3) + t1 * qJ(4) - t2 * t6 + t31 * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t17 ^ 2, -t3 * t17 + t4;];
tauc_reg  = t9;
