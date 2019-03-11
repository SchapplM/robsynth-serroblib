% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% tauc_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:24
% EndTime: 2019-03-08 18:12:24
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->15), mult. (49->19), div. (0->0), fcn. (21->2), ass. (0->13)
t4 = sin(qJ(3));
t9 = qJD(2) * t4;
t3 = qJD(3) * qJ(4) + t9;
t5 = cos(qJ(3));
t12 = t3 * t5;
t6 = qJD(3) ^ 2;
t11 = t6 * t4;
t10 = t6 * t5;
t8 = qJD(2) * t5;
t7 = qJD(3) * pkin(3);
t2 = qJD(4) - t7 - t8;
t1 = (qJD(4) + t8) * qJD(3);
t13 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t11, -t10, -t11, t10, t1 * t4 + (t12 + (t2 - t8) * t4) * qJD(3); 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * qJD(4), t1 * qJ(4) + t3 * qJD(4) + (-t12 + (-t2 - t7) * t4) * qJD(2); 0, 0, 0, 0, 0, 0, -t6 (-t3 + t9) * qJD(3);];
tauc_reg  = t13;
