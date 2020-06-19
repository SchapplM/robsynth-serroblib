% Calculate vector of inverse dynamics joint torques for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
% MDP [3x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S1R1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S1R1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_invdynJ_fixb_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1R1_invdynJ_fixb_mdp_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1R1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1R1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_invdynJ_fixb_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [3 1]), ...
  'S1R1_invdynJ_fixb_mdp_slag_vp: MDP has to be [3x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:57
% EndTime: 2020-06-19 09:12:57
% DurationCPUTime: 0.03s
% Computational Cost: add. (4->4), mult. (7->7), div. (0->0), fcn. (4->2), ass. (0->3)
t4 = cos(qJ(1));
t3 = sin(qJ(1));
t1 = [qJDD(1) * MDP(1) + (g(1) * t3 - g(2) * t4) * MDP(2) + (g(1) * t4 + g(2) * t3) * MDP(3);];
tau = t1;
