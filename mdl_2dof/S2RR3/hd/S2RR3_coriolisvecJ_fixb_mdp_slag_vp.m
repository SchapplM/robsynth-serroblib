% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S2RR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S2RR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:28
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.05s
% Computational Cost: add. (10->4), mult. (28->7), div. (0->0), fcn. (8->2), ass. (0->3)
t13 = pkin(1) * (cos(qJ(2)) * MDP(6) + sin(qJ(2)) * MDP(5));
t6 = qJD(1) + qJD(2);
t1 = [(-qJD(1) - t6) * qJD(2) * t13; (-qJD(2) + t6) * qJD(1) * t13;];
tauc = t1;
