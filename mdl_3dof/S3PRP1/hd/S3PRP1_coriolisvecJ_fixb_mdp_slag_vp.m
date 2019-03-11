% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3PRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:01
% EndTime: 2019-03-08 18:03:01
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->17), mult. (58->25), div. (0->0), fcn. (21->2), ass. (0->11)
t22 = qJD(2) * pkin(2);
t16 = sin(qJ(2));
t21 = t16 * qJD(1);
t17 = cos(qJ(2));
t20 = t17 * qJD(1);
t19 = MDP(7) * qJD(2);
t18 = qJD(2) ^ 2;
t15 = qJD(2) * qJ(3) + t21;
t14 = qJD(3) - t20 - t22;
t13 = (qJD(3) + t20) * qJD(2);
t1 = [(t15 * t19 + (-MDP(4) + MDP(6)) * t18) * t17 + ((t13 + (t14 - t20) * qJD(2)) * MDP(7) + (-MDP(3) - MDP(5)) * t18) * t16; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (t13 * qJ(3) + t15 * qJD(3) + (-t15 * t17 + (-t14 - t22) * t16) * qJD(1)) * MDP(7); -t18 * MDP(6) + (-t15 + t21) * t19;];
tauc  = t1;
