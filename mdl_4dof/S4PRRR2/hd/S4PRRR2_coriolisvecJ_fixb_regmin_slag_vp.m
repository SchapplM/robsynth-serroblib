% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% tauc_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:24
% EndTime: 2019-07-18 13:27:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (78->30), mult. (188->50), div. (0->0), fcn. (90->4), ass. (0->26)
t21 = qJD(3) + qJD(4);
t7 = sin(qJ(4));
t8 = sin(qJ(3));
t28 = t7 * t8;
t9 = cos(qJ(4));
t27 = t8 * t9;
t25 = pkin(1) * qJD(2);
t18 = t25 * t28;
t26 = t21 * t18;
t24 = qJD(4) * t8;
t10 = cos(qJ(3));
t23 = qJD(3) * t10;
t22 = t10 * qJD(2);
t20 = pkin(1) * t22;
t6 = qJD(2) + qJD(3);
t1 = t6 * pkin(2) + t20;
t5 = qJD(2) + t21;
t19 = (-qJD(4) + t5) * t1;
t17 = -t10 * t7 - t27;
t16 = (-qJD(3) + t6) * t25;
t15 = pkin(1) * qJD(3) * (-qJD(2) - t6);
t14 = qJD(4) * (-pkin(2) * t5 - t1);
t13 = qJD(4) * (-(t10 * pkin(1) + pkin(2)) * t5 - t1);
t12 = (t10 * t9 - t28) * t5;
t11 = t17 * qJD(3) - t9 * t24;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t8 * t15, t10 * t15, 0, t7 * t13 + (qJD(2) + t5) * pkin(1) * t11, t9 * t13 + (t7 * t5 * t24 + (-t9 * t22 - t12) * qJD(3)) * pkin(1) + t26; 0, 0, 0, 0, 0, t8 * t16, t10 * t16, 0, t7 * t14 + (-t17 * t5 + t11) * t25, t9 * t14 + (-t9 * t23 + t12) * t25 + t26; 0, 0, 0, 0, 0, 0, 0, 0, t7 * t19 + (-t7 * t23 + (t5 - t21) * t27) * t25, -t5 * t18 + (-qJD(3) * t20 + t19) * t9 + t26;];
tauc_reg  = t2;
