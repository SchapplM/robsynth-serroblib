% Calculate vector of inverse dynamics joint torques for
% S4PPRP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:40
% EndTime: 2019-03-08 18:14:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (73->38), mult. (119->49), div. (0->0), fcn. (70->2), ass. (0->14)
t29 = sin(qJ(3));
t30 = cos(qJ(3));
t36 = -t29 * qJD(1) + t30 * qJD(2);
t34 = qJDD(1) - g(2);
t33 = qJDD(2) - g(1);
t25 = t30 * qJD(1) + t29 * qJD(2);
t31 = qJD(3) ^ 2;
t26 = t30 * qJDD(2);
t23 = -t29 * qJDD(3) - t30 * t31;
t22 = -t30 * qJDD(3) + t29 * t31;
t21 = qJD(3) * pkin(3) + t36;
t20 = qJD(3) * t36 + t30 * qJDD(1) + t29 * qJDD(2);
t19 = qJDD(3) * pkin(3) - t25 * qJD(3) - t29 * qJDD(1) + t26;
t1 = [t23 * MDP(4) + t22 * MDP(5) + (-t19 * t29 + t20 * t30 - g(2) + (-t21 * t30 - t25 * t29) * qJD(3)) * MDP(6) + (MDP(1) + MDP(2)) * t34; t33 * MDP(2) - t22 * MDP(4) + t23 * MDP(5) + (t19 * t30 + t20 * t29 - g(1) + (-t21 * t29 + t25 * t30) * qJD(3)) * MDP(6); qJDD(3) * MDP(3) + t26 * MDP(4) + (-g(1) * MDP(4) - t34 * MDP(5)) * t30 + (-t34 * MDP(4) - t33 * MDP(5)) * t29 + ((t21 - t36) * t25 + (-g(1) * t30 + g(2) * t29 + t19) * pkin(3)) * MDP(6); (qJDD(4) + g(3)) * MDP(6);];
tau  = t1;
