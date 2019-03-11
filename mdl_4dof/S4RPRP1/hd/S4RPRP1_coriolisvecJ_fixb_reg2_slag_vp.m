% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:43
% EndTime: 2019-03-08 18:29:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (147->32), mult. (361->42), div. (0->0), fcn. (180->4), ass. (0->28)
t16 = cos(pkin(6)) * pkin(1) + pkin(2);
t21 = sin(qJ(3));
t22 = cos(qJ(3));
t35 = sin(pkin(6)) * pkin(1);
t33 = t21 * t16 + t22 * t35;
t13 = t16 * qJD(1);
t29 = qJD(1) * t35;
t6 = t21 * t13 + t22 * t29;
t36 = t6 * qJD(3);
t28 = t21 * t29;
t31 = qJD(3) * t22;
t34 = qJD(3) * t28 - t13 * t31;
t5 = t22 * t13 - t28;
t32 = qJD(4) - t5;
t30 = t21 * t35;
t18 = qJD(1) + qJD(3);
t17 = t18 * qJD(4);
t1 = t17 - t34;
t27 = t5 * t18 + t34;
t8 = -qJD(3) * t30 + t16 * t31;
t26 = t22 * t16 - t30;
t24 = t6 * t18 - t36;
t9 = t33 * qJD(3);
t23 = -t9 * t18 - t36;
t7 = qJD(4) + t8;
t3 = t18 * qJ(4) + t6;
t2 = -t18 * pkin(3) + t32;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t8 * t18 + t34, 0, -t26 * t36 - t34 * t33 - t5 * t9 + t6 * t8, 0, 0, 0, 0, 0, 0, t23, 0, t7 * t18 + t1, t1 * (qJ(4) + t33) + t3 * t7 + t36 * (-pkin(3) - t26) + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t27, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0.2e1 * t17 - t27, -pkin(3) * t36 + t1 * qJ(4) - t2 * t6 + t32 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18 ^ 2, -t3 * t18 + t36;];
tauc_reg  = t4;
