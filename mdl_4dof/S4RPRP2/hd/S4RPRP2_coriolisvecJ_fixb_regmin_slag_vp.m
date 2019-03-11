% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauc_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:53
% EndTime: 2019-03-08 18:30:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (123->30), mult. (206->43), div. (0->0), fcn. (77->2), ass. (0->25)
t13 = sin(qJ(3));
t27 = t13 * qJD(2);
t14 = cos(qJ(3));
t28 = qJD(3) * t14;
t15 = -pkin(1) - pkin(2);
t10 = t15 * qJD(1) + qJD(2);
t29 = t13 * t10;
t2 = -(qJ(2) * t28 + t27) * qJD(1) - qJD(3) * t29;
t24 = qJD(1) - qJD(3);
t26 = qJ(2) * qJD(1);
t7 = t14 * t26 + t29;
t30 = -t7 * t24 + t2;
t25 = qJD(1) * qJD(2);
t22 = 0.2e1 * t25;
t21 = t13 * t26;
t20 = t24 ^ 2;
t6 = t14 * t10 - t21;
t19 = t14 * qJ(2) + t13 * t15;
t18 = -t13 * qJ(2) + t14 * t15;
t1 = -qJD(3) * t21 + t10 * t28 + t14 * t25;
t16 = qJD(1) ^ 2;
t5 = -t19 * qJD(3) - t27;
t4 = t14 * qJD(2) + t18 * qJD(3);
t3 = -pkin(3) * t24 + t6;
t8 = [0, 0, 0, 0, t22, qJ(2) * t22, 0, -t24 * t5 - t2, t24 * t4 + t1, t1 * t19 + t7 * t4 + t2 * (-pkin(3) + t18) + t3 * t5; 0, 0, 0, 0, -t16, -t16 * qJ(2), 0, -t13 * t20, -t14 * t20, t30 * t14 + (t24 * t3 + t1) * t13; 0, 0, 0, 0, 0, 0, 0, t30, -t24 * t6 - t1, t2 * pkin(3) + (t3 - t6) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t8;
