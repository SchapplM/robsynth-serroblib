% Calculate inertial parameters regressor of coriolis joint torque vector for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% tauc_reg [2x(2*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S2RR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:02
% EndTime: 2019-03-08 18:00:03
% DurationCPUTime: 0.07s
% Computational Cost: add. (5->4), mult. (38->17), div. (0->0), fcn. (16->2), ass. (0->12)
t3 = sin(qJ(2));
t4 = cos(qJ(2));
t13 = t3 * t4;
t5 = qJD(2) ^ 2;
t12 = t5 * t3;
t11 = t5 * t4;
t10 = t3 ^ 2 - t4 ^ 2;
t6 = qJD(1) ^ 2;
t9 = t6 * t13;
t8 = qJD(1) * qJD(2);
t7 = t8 * t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t7, -0.2e1 * t10 * t8, -t11, -0.2e1 * t7, t12, 0, pkin(1) * t11, -pkin(1) * t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t10 * t6, 0, t9, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t1;
