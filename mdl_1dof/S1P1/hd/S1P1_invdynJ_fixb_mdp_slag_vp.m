% Calculate vector of inverse dynamics joint torques for
% S1P1
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
%   pkin=[dummy]';
% MDP [1x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S1P1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S1P1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [1 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vp: MDP has to be [1x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:56
% EndTime: 2021-01-14 12:21:56
% DurationCPUTime: 0.02s
% Computational Cost: add. (1->1), mult. (1->1), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [(qJDD(1) - g(3)) * MDP(1);];
tau = t1;
