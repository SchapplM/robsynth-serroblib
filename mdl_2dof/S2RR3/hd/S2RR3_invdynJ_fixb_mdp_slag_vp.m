% Calculate vector of inverse dynamics joint torques for
% S2RR3
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (43->22), mult. (57->29), div. (0->0), fcn. (26->6), ass. (0->13)
t28 = sin(qJ(2));
t30 = cos(qJ(2));
t38 = t28 * MDP(5) + t30 * MDP(6);
t36 = t28 * MDP(6);
t27 = qJ(1) + qJ(2);
t23 = sin(t27);
t24 = cos(t27);
t25 = qJDD(1) + qJDD(2);
t32 = t25 * MDP(4) + (t30 * qJDD(1) * pkin(1) + g(1) * t23 - g(2) * t24) * MDP(5) + (g(1) * t24 + g(2) * t23) * MDP(6);
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t26 = qJD(1) + qJD(2);
t1 = [qJDD(1) * MDP(1) + (g(1) * t29 - g(2) * t31) * MDP(2) + (g(1) * t31 + g(2) * t29) * MDP(3) + (t25 * t30 * MDP(5) + (-qJDD(1) - t25) * t36 + t38 * qJD(2) * (-qJD(1) - t26)) * pkin(1) + t32; (-qJDD(1) * t36 + t38 * qJD(1) * (-qJD(2) + t26)) * pkin(1) + t32;];
tau = t1;
