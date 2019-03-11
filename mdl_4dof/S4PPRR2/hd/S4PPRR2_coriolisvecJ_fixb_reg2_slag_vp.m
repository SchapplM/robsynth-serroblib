% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:01
% DurationCPUTime: 0.17s
% Computational Cost: add. (207->42), mult. (588->63), div. (0->0), fcn. (498->6), ass. (0->32)
t22 = sin(pkin(6));
t23 = cos(pkin(6));
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t39 = -t25 * t22 + t27 * t23;
t14 = t39 * qJD(1);
t19 = t27 * t22 + t25 * t23;
t15 = t19 * qJD(1);
t11 = qJD(3) * pkin(3) + t14;
t16 = t39 * qJD(3);
t12 = qJD(1) * t16;
t26 = cos(qJ(4));
t17 = t19 * qJD(3);
t13 = qJD(1) * t17;
t24 = sin(qJ(4));
t38 = t24 * t15;
t32 = -qJD(4) * t38 - t24 * t13;
t1 = (qJD(4) * t11 + t12) * t26 + t32;
t36 = t26 * t15;
t21 = qJD(3) + qJD(4);
t34 = -pkin(3) * t21 - t11;
t33 = -t24 * t12 - t26 * t13;
t6 = t24 * t11 + t36;
t30 = -t24 * t19 + t26 * t39;
t29 = t26 * t19 + t24 * t39;
t2 = -qJD(4) * t6 + t33;
t8 = t26 * t14 - t38;
t7 = -t24 * t14 - t36;
t5 = t26 * t11 - t38;
t4 = -t29 * qJD(4) - t24 * t16 - t26 * t17;
t3 = t30 * qJD(4) + t26 * t16 - t24 * t17;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * qJD(3), -t16 * qJD(3), 0, t12 * t19 - t13 * t39 - t14 * t17 + t15 * t16, 0, 0, 0, 0, 0, 0, t4 * t21, -t3 * t21, 0, t1 * t29 + t2 * t30 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t21 + (t34 * t24 - t36) * qJD(4) + t33, t8 * t21 + (t34 * qJD(4) - t12) * t26 - t32, 0, -t5 * t7 - t6 * t8 + (t1 * t24 + t2 * t26 + (-t24 * t5 + t26 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t21 + t2, t5 * t21 - t1, 0, 0;];
tauc_reg  = t9;
