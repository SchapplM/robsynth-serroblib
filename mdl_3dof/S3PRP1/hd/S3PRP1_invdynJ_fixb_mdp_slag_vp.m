% Calculate vector of inverse dynamics joint torques for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3PRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:01
% EndTime: 2019-03-08 18:03:01
% DurationCPUTime: 0.09s
% Computational Cost: add. (67->45), mult. (102->52), div. (0->0), fcn. (48->2), ass. (0->17)
t33 = sin(qJ(2));
t34 = cos(qJ(2));
t44 = g(1) * t33 - g(2) * t34;
t42 = qJDD(2) * pkin(2);
t41 = t33 * qJD(1);
t40 = t34 * qJD(1);
t39 = qJDD(2) * qJ(3);
t31 = t34 * qJDD(1);
t38 = t31 + t44;
t28 = qJD(2) * qJ(3) + t41;
t37 = (-qJD(2) * pkin(2) + qJD(3) - t40) * t33 + t28 * t34;
t30 = t33 * qJDD(1);
t36 = g(1) * t34 + g(2) * t33 - t30;
t24 = qJD(2) * t41 + qJDD(3) - t31 - t42;
t35 = qJD(2) ^ 2;
t23 = t39 + t30 + (qJD(3) + t40) * qJD(2);
t1 = [(qJDD(1) - g(2)) * MDP(1) + (qJD(2) * t37 + t23 * t33 - t24 * t34 - g(2)) * MDP(7) + (MDP(3) + MDP(5)) * (qJDD(2) * t34 - t33 * t35) + (-MDP(4) + MDP(6)) * (qJDD(2) * t33 + t34 * t35); qJDD(2) * MDP(2) + t38 * MDP(3) + t36 * MDP(4) + (-qJDD(3) + t38 + 0.2e1 * t42) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) - t36 + 0.2e1 * t39) * MDP(6) + (t23 * qJ(3) + t28 * qJD(3) - t24 * pkin(2) - g(1) * (-pkin(2) * t33 + qJ(3) * t34) - g(2) * (pkin(2) * t34 + qJ(3) * t33) - t37 * qJD(1)) * MDP(7); -qJDD(2) * MDP(5) - t35 * MDP(6) + (-qJD(2) * t28 + t24 - t44) * MDP(7);];
tau  = t1;
