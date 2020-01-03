% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPR4
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:01
% DurationCPUTime: 0.13s
% Computational Cost: add. (118->27), mult. (243->50), div. (0->0), fcn. (98->2), ass. (0->25)
t19 = 2 * qJD(2);
t10 = cos(qJ(4));
t11 = (-pkin(2) - pkin(5));
t5 = t11 * qJD(2) + qJD(3);
t9 = sin(qJ(4));
t3 = -t9 * qJD(1) + t10 * t5;
t1 = t3 * qJD(4);
t4 = t10 * qJD(1) + t9 * t5;
t2 = t4 * qJD(4);
t14 = (t10 * t4 - t3 * t9) * qJD(4) + t1 * t9 - t2 * t10;
t28 = -t10 ^ 2 + t9 ^ 2;
t27 = t10 * t9;
t12 = qJD(4) ^ 2;
t26 = t12 * t9;
t25 = t12 * t10;
t13 = qJD(2) ^ 2;
t24 = -t12 - t13;
t23 = t13 * qJ(3);
t22 = qJ(3) * qJD(4);
t21 = qJD(2) * qJD(4);
t20 = t13 * t27;
t18 = qJD(3) * t19;
t17 = t21 * t27;
t6 = qJ(3) * t18;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t26, 0, t1 * t10 + t2 * t9 + (-t10 * t3 - t4 * t9) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t6, -0.2e1 * t17, 0.2e1 * t28 * t21, -t26, 0.2e1 * t17, -t25, 0, -t11 * t26 + (qJD(3) * t9 + t10 * t22) * t19, -t11 * t25 + (qJD(3) * t10 - t9 * t22) * t19, -t14, t14 * t11 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t23, 0, 0, 0, 0, 0, 0, t24 * t9, t24 * t10, 0, t14 - t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t28 * t13, 0, -t20, 0, 0, -t10 * t23, t9 * t23, 0, 0;];
tauc_reg = t7;
