% Calculate vector of inverse dynamics joint torques for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:11
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->26), mult. (86->46), div. (0->0), fcn. (45->4), ass. (0->15)
t21 = sin(qJ(2));
t23 = cos(qJ(2));
t33 = t21 * t23;
t19 = t21 ^ 2;
t32 = -t23 ^ 2 + t19;
t31 = qJD(1) * qJD(2);
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t30 = g(1) * t22 + g(3) * t24;
t25 = qJD(2) ^ 2;
t29 = qJDD(2) * t23 - t25 * t21;
t28 = qJDD(2) * t21 + t25 * t23;
t27 = pkin(1) * qJDD(1) + t30;
t26 = qJD(1) ^ 2;
t1 = [qJDD(1) * MDP(1) + t30 * MDP(3) + (t19 * qJDD(1) + 0.2e1 * t31 * t33) * MDP(4) + 0.2e1 * (qJDD(1) * t33 - t32 * t31) * MDP(5) - t28 * MDP(6) - t29 * MDP(7) + (t29 * MDP(10) + t28 * MDP(9)) * pkin(1) + (-t21 * MDP(10) + t23 * MDP(9) + MDP(2)) * (-g(1) * t24 + g(3) * t22); qJDD(2) * MDP(8) + t32 * MDP(5) * t26 + (t27 * MDP(10) - qJDD(1) * MDP(7) + g(2) * MDP(9)) * t23 + (-t26 * t23 * MDP(4) - g(2) * MDP(10) - qJDD(1) * MDP(6) + t27 * MDP(9)) * t21;];
tau = t1;
