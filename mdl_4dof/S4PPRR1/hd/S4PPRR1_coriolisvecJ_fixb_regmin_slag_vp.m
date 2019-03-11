% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:49
% EndTime: 2019-03-08 18:15:49
% DurationCPUTime: 0.12s
% Computational Cost: add. (45->18), mult. (102->30), div. (0->0), fcn. (68->4), ass. (0->18)
t4 = qJD(3) + qJD(4);
t5 = sin(qJ(4));
t6 = sin(qJ(3));
t21 = t5 * t6;
t8 = cos(qJ(3));
t20 = t5 * t8;
t7 = cos(qJ(4));
t19 = t6 * t7;
t18 = t7 * t8;
t16 = qJD(2) * qJD(3);
t3 = qJD(3) * pkin(3) + qJD(2) * t8;
t14 = (-qJD(4) + t4) * t3;
t13 = -t19 - t20;
t12 = t18 - t21;
t11 = qJD(4) * (-pkin(3) * t4 - t3);
t10 = t13 * qJD(3);
t9 = qJD(3) ^ 2;
t1 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t9 * t6, -t9 * t8, 0 (qJD(4) * t13 + t10) * t4, -t4 ^ 2 * t12; 0, 0, 0, 0, 0, 0, t5 * t11 + (-qJD(4) * t19 - t13 * t4 + t10) * qJD(2), t7 * t11 + (-qJD(3) * t18 + (t12 + t21) * t4) * qJD(2); 0, 0, 0, 0, 0, 0, t14 * t5 - t16 * t20 (-t16 * t8 + t14) * t7;];
tauc_reg  = t1;
