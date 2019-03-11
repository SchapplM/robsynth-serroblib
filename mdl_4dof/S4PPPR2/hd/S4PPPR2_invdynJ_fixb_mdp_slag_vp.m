% Calculate vector of inverse dynamics joint torques for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:18
% EndTime: 2019-03-08 18:10:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (52->29), mult. (88->42), div. (0->0), fcn. (68->4), ass. (0->12)
t28 = qJD(4) ^ 2;
t24 = sin(pkin(5));
t32 = qJDD(1) * t24 ^ 2 - g(2);
t31 = qJDD(1) - g(2);
t30 = t24 * qJDD(1);
t25 = cos(pkin(5));
t26 = sin(qJ(4));
t27 = cos(qJ(4));
t18 = -t24 * t27 + t25 * t26;
t29 = t24 * t26 + t25 * t27;
t19 = -qJDD(1) * t25 + qJDD(3);
t1 = [t31 * MDP(1) + (qJDD(1) * t25 ^ 2 + t32) * MDP(2) + (-t19 * t25 + t32) * MDP(3) + (-qJDD(4) * t29 + t18 * t28) * MDP(5) + (t18 * qJDD(4) + t28 * t29) * MDP(6); (MDP(2) + MDP(3)) * (qJDD(2) - g(3)); (-g(1) * t24 - t25 * t31 + qJDD(3)) * MDP(3) + (qJDD(4) * t27 - t26 * t28) * MDP(5) + (-qJDD(4) * t26 - t27 * t28) * MDP(6); qJDD(4) * MDP(4) + (g(1) * t18 + g(2) * t29 + t27 * t19 - t26 * t30) * MDP(5) + (g(1) * t29 - g(2) * t18 - t26 * t19 - t27 * t30) * MDP(6);];
tau  = t1;
