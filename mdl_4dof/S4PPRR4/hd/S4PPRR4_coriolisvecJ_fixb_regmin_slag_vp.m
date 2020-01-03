% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [4x12]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (73->21), mult. (231->44), div. (0->0), fcn. (164->6), ass. (0->26)
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t10 = t18 * t13 + t16 * t14;
t8 = t10 * qJD(3);
t9 = t16 * t13 - t18 * t14;
t33 = t9 * qJD(1);
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t32 = t15 * t17;
t19 = qJD(4) ^ 2;
t31 = t19 * t15;
t30 = t19 * t17;
t29 = t15 ^ 2 - t17 ^ 2;
t28 = qJD(3) * pkin(3);
t27 = t15 * qJD(4);
t26 = t17 * qJD(4);
t25 = 0.2e1 * qJD(3) * qJD(4);
t1 = t33 - t28;
t7 = t9 * qJD(3);
t24 = qJD(1) * t7 - t1 * qJD(3);
t23 = pkin(5) * t19;
t22 = qJD(4) * (t1 - t33 - t28);
t20 = qJD(3) ^ 2;
t2 = [0, 0, 0, -t8 * qJD(3), t7 * qJD(3), 0, 0, 0, 0, 0, t7 * t27 - t10 * t30 + (-t17 * t8 + t9 * t27) * qJD(3), t7 * t26 + t10 * t31 + (t15 * t8 + t9 * t26) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t30; 0, 0, 0, 0, 0, t25 * t32, -t29 * t25, t30, -t31, 0, t15 * t22 - t23 * t17, t23 * t15 + t17 * t22; 0, 0, 0, 0, 0, -t20 * t32, t29 * t20, 0, 0, 0, t24 * t15, t24 * t17;];
tauc_reg = t2;
