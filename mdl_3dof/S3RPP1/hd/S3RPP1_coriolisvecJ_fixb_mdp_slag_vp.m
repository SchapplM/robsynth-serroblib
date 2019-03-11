% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3RPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:02
% EndTime: 2019-03-08 18:05:02
% DurationCPUTime: 0.03s
% Computational Cost: add. (20->14), mult. (41->19), div. (0->0), fcn. (0->0), ass. (0->7)
t15 = (qJ(2) * MDP(6) + MDP(5) + MDP(7));
t13 = (-pkin(1) - qJ(3));
t11 = (MDP(9) * qJD(1));
t9 = (qJD(1) ^ 2);
t8 = (qJD(1) * qJ(2) + qJD(3));
t7 = (t13 * qJD(1) + qJD(2));
t1 = [(t8 * qJD(2) - t7 * qJD(3) + (qJ(2) * qJD(2) - t13 * qJD(3)) * qJD(1)) * MDP(9) + 2 * (qJD(3) * MDP(8) + t15 * qJD(2)) * qJD(1); (-qJD(3) - t8) * t11 - t15 * t9; -t9 * MDP(8) + (qJD(2) + t7) * t11;];
tauc  = t1;
