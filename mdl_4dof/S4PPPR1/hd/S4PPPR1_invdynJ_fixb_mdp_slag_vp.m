% Calculate vector of inverse dynamics joint torques for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:09:05
% EndTime: 2019-03-08 18:09:05
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->23), mult. (61->31), div. (0->0), fcn. (42->4), ass. (0->11)
t23 = MDP(2) + MDP(3);
t22 = qJD(4) ^ 2;
t21 = cos(qJ(4));
t20 = sin(qJ(4));
t19 = cos(pkin(5));
t18 = sin(pkin(5));
t16 = -t20 * qJDD(4) - t21 * t22;
t15 = -t21 * qJDD(4) + t20 * t22;
t13 = t18 * t21 + t19 * t20;
t12 = -t18 * t20 + t19 * t21;
t1 = [(MDP(1) + t23) * (qJDD(1) - g(3)); t16 * MDP(5) + t15 * MDP(6) + t23 * (-g(1) * t18 + g(2) * t19 + qJDD(2)); (-g(1) * t19 - g(2) * t18 + qJDD(3)) * MDP(3) - t15 * MDP(5) + t16 * MDP(6); qJDD(4) * MDP(4) + (-g(1) * t12 - g(2) * t13 - t20 * qJDD(2) + t21 * qJDD(3)) * MDP(5) + (g(1) * t13 - g(2) * t12 - t21 * qJDD(2) - t20 * qJDD(3)) * MDP(6);];
tau  = t1;
