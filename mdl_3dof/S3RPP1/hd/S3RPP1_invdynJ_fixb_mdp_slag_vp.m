% Calculate vector of inverse dynamics joint torques for
% S3RPP1
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
%   pkin=[a2,a3,d1]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:02
% EndTime: 2019-03-08 18:05:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (77->51), mult. (100->54), div. (0->0), fcn. (28->2), ass. (0->21)
t38 = -pkin(1) - qJ(3);
t51 = (qJDD(1) * t38);
t36 = (qJD(1) * qJD(2));
t31 = 2 * t36;
t39 = sin(qJ(1));
t40 = cos(qJ(1));
t50 = t40 * pkin(1) + t39 * qJ(2);
t49 = g(1) * t39 - g(2) * t40;
t48 = qJDD(1) * pkin(1);
t37 = qJ(2) * qJDD(1);
t47 = qJD(3) * qJD(1);
t46 = qJDD(2) - t49;
t45 = qJDD(3) + t36 + t37;
t44 = qJDD(2) - t47;
t43 = g(1) * t40 + g(2) * t39;
t42 = t31 + 0.2e1 * t37 - t43;
t41 = qJD(1) ^ 2;
t30 = t40 * qJ(2);
t28 = qJ(2) * qJD(1) + qJD(3);
t27 = t38 * qJD(1) + qJD(2);
t1 = [qJDD(1) * MDP(1) + t49 * MDP(2) + t43 * MDP(3) + (t46 - 0.2e1 * t48) * MDP(4) + t42 * MDP(5) + (-(qJDD(2) - t48) * pkin(1) - g(1) * (-pkin(1) * t39 + t30) - g(2) * t50 + (t31 + t37) * qJ(2)) * MDP(6) + (qJDD(3) + t42) * MDP(7) + (-t46 + 0.2e1 * t47 - (2 * t51)) * MDP(8) + ((t44 + t51) * t38 - t27 * qJD(3) + t45 * qJ(2) + t28 * qJD(2) - g(1) * (t38 * t39 + t30) - g(2) * (qJ(3) * t40 + t50)) * MDP(9); t46 * MDP(6) + (-qJD(1) * t28 + t44 - t49) * MDP(9) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t41 + (-pkin(1) * MDP(6) + t38 * MDP(9) + MDP(4) - MDP(8)) * qJDD(1); qJDD(1) * MDP(7) - t41 * MDP(8) + (qJD(1) * t27 - t43 + t45) * MDP(9);];
tau  = t1;
