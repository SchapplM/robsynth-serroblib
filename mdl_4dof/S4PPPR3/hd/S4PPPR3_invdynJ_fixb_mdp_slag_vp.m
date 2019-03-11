% Calculate vector of inverse dynamics joint torques for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:26
% EndTime: 2019-03-08 18:11:26
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->29), mult. (103->37), div. (0->0), fcn. (84->6), ass. (0->15)
t36 = qJD(4) ^ 2;
t30 = sin(pkin(5));
t31 = cos(pkin(5));
t32 = sin(qJ(4));
t33 = cos(qJ(4));
t21 = -t33 * t30 - t32 * t31;
t35 = -t32 * t30 + t33 * t31;
t34 = t21 * qJDD(4) - t35 * t36;
t28 = pkin(5) + qJ(4);
t27 = cos(t28);
t26 = sin(t28);
t23 = t31 * qJDD(1) + t30 * qJDD(2);
t22 = -t30 * qJDD(1) + t31 * qJDD(2);
t19 = t35 * qJDD(4) + t21 * t36;
t1 = [(-t22 * t30 + t23 * t31 - g(2)) * MDP(3) + t34 * MDP(5) - t19 * MDP(6) + (MDP(1) + MDP(2)) * (qJDD(1) - g(2)); (qJDD(2) - g(1)) * MDP(2) + (t22 * t31 + t23 * t30 - g(1)) * MDP(3) + t19 * MDP(5) + t34 * MDP(6); (qJDD(3) + g(3)) * MDP(3); qJDD(4) * MDP(4) + (-g(1) * t27 + g(2) * t26 + t33 * t22 - t32 * t23) * MDP(5) + (g(1) * t26 + g(2) * t27 - t32 * t22 - t33 * t23) * MDP(6);];
tau  = t1;
