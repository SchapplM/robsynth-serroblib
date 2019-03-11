% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:56
% EndTime: 2019-03-08 18:17:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (20->14), mult. (41->19), div. (0->0), fcn. (0->0), ass. (0->7)
t15 = (qJ(3) * MDP(7) + MDP(6) + MDP(8));
t13 = (-pkin(2) - qJ(4));
t11 = (MDP(10) * qJD(2));
t9 = (qJD(2) ^ 2);
t8 = (qJD(2) * qJ(3) + qJD(4));
t7 = (t13 * qJD(2) + qJD(3));
t1 = [0; (t8 * qJD(3) - t7 * qJD(4) + (qJ(3) * qJD(3) - t13 * qJD(4)) * qJD(2)) * MDP(10) + 2 * (qJD(4) * MDP(9) + t15 * qJD(3)) * qJD(2); (-qJD(4) - t8) * t11 - t15 * t9; -t9 * MDP(9) + (qJD(3) + t7) * t11;];
tauc  = t1;
