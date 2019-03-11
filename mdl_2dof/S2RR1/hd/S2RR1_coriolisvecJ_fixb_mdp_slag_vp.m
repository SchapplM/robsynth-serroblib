% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S2RR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S2RR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:12
% EndTime: 2019-03-08 18:00:12
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->6), mult. (39->16), div. (0->0), fcn. (12->2), ass. (0->4)
t14 = sin(qJ(2));
t15 = cos(qJ(2));
t20 = -t14 * t15 * MDP(4) + (t14 ^ 2 - t15 ^ 2) * MDP(5);
t1 = [(-0.2e1 * t20 * qJD(1) + (-t15 * MDP(6) + t14 * MDP(7) + (-MDP(10) * t14 + MDP(9) * t15) * pkin(1)) * qJD(2)) * qJD(2); t20 * qJD(1) ^ 2;];
tauc  = t1;
