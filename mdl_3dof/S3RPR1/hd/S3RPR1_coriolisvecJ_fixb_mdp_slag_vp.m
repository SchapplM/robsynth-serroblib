% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3RPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (49->22), mult. (86->39), div. (0->0), fcn. (26->2), ass. (0->13)
t18 = sin(qJ(3));
t19 = cos(qJ(3));
t29 = t19 * MDP(9);
t32 = t18 * MDP(8) + t29;
t31 = (qJ(2) * MDP(6) + MDP(5));
t20 = -pkin(1) - pkin(2);
t25 = qJD(1) - qJD(3);
t27 = qJD(3) + t25;
t26 = qJD(1) * qJ(2);
t23 = t27 * MDP(8);
t22 = qJD(2) * (qJD(1) + t25);
t16 = t20 * qJD(1) + qJD(2);
t1 = [(t18 * t22 + (-(-qJ(2) * t19 - t18 * t20) * t25 + t19 * t26 + t18 * t16) * qJD(3)) * MDP(8) + (t19 * t22 + ((-qJ(2) * t18 + t19 * t20) * t25 - t18 * t26 + t19 * t16) * qJD(3)) * MDP(9) + (2 * t31 * qJD(1) * qJD(2)); -(t31 * qJD(1) ^ 2) - t25 ^ 2 * t32; (-t18 * t23 - t27 * t29) * t16 + (-t32 * qJD(2) + (t27 * MDP(9) * t18 - t19 * t23) * qJ(2)) * qJD(1);];
tauc  = t1;
