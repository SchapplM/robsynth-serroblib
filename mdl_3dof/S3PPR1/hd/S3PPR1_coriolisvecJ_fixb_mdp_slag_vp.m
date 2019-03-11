% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% MDP [5x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3PPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [5 1]), ...
  'S3PPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [5x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:02:01
% EndTime: 2019-03-08 18:02:01
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->2), mult. (28->7), div. (0->0), fcn. (12->2), ass. (0->4)
t8 = qJD(3) ^ 2;
t7 = cos(qJ(3));
t6 = sin(qJ(3));
t1 = [(-MDP(4) * t7 + MDP(5) * t6) * t8; (-MDP(4) * t6 - MDP(5) * t7) * t8; 0;];
tauc  = t1;
