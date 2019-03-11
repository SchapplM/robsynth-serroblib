% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% tauc_reg [4x6]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:50
% EndTime: 2019-03-08 18:14:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->14), mult. (76->22), div. (0->0), fcn. (46->2), ass. (0->11)
t7 = sin(qJ(3));
t8 = cos(qJ(3));
t14 = -qJD(1) * t7 + t8 * qJD(2);
t9 = qJD(3) ^ 2;
t13 = t7 * t9;
t12 = t8 * t9;
t5 = qJD(1) * t8 + qJD(2) * t7;
t3 = qJD(3) * t5;
t2 = t14 * qJD(3);
t1 = qJD(3) * pkin(3) + t14;
t4 = [0, 0, 0, -t12, t13, t2 * t8 + t3 * t7 + (-t1 * t8 - t5 * t7) * qJD(3); 0, 0, 0, -t13, -t12, t2 * t7 - t3 * t8 + (-t1 * t7 + t5 * t8) * qJD(3); 0, 0, 0, 0, 0, -pkin(3) * t3 + (t1 - t14) * t5; 0, 0, 0, 0, 0, 0;];
tauc_reg  = t4;
