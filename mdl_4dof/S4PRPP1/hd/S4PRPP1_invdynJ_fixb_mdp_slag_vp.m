% Calculate vector of inverse dynamics joint torques for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:06
% EndTime: 2019-03-08 18:18:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (110->55), mult. (103->55), div. (0->0), fcn. (28->2), ass. (0->22)
t44 = -pkin(2) - qJ(4);
t55 = (t44 * qJDD(2));
t41 = (qJD(2) * qJD(3));
t38 = 2 * t41;
t40 = pkin(5) + qJ(2);
t36 = sin(t40);
t37 = cos(t40);
t54 = t37 * pkin(2) + t36 * qJ(3);
t53 = g(1) * t36 - g(2) * t37;
t52 = qJDD(2) * pkin(2);
t42 = qJ(3) * qJDD(2);
t51 = qJD(4) * qJD(2);
t50 = qJDD(3) - t53;
t49 = qJDD(4) + t41 + t42;
t48 = qJDD(3) - t51;
t47 = g(1) * t37 + g(2) * t36;
t46 = t38 + 0.2e1 * t42 - t47;
t45 = qJD(2) ^ 2;
t35 = qJ(3) * qJD(2) + qJD(4);
t31 = t37 * qJ(3);
t29 = t44 * qJD(2) + qJD(3);
t1 = [(MDP(1) + MDP(7) + MDP(10)) * (qJDD(1) - g(3)); qJDD(2) * MDP(2) + t53 * MDP(3) + t47 * MDP(4) + (t50 - 0.2e1 * t52) * MDP(5) + t46 * MDP(6) + (-(qJDD(3) - t52) * pkin(2) - g(1) * (-t36 * pkin(2) + t31) - g(2) * t54 + (t38 + t42) * qJ(3)) * MDP(7) + (qJDD(4) + t46) * MDP(8) + (-t50 + 0.2e1 * t51 - (2 * t55)) * MDP(9) + ((t48 + t55) * t44 - t29 * qJD(4) + t49 * qJ(3) + t35 * qJD(3) - g(1) * (t44 * t36 + t31) - g(2) * (t37 * qJ(4) + t54)) * MDP(10); t50 * MDP(7) + (-t35 * qJD(2) + t48 - t53) * MDP(10) + (-qJ(3) * MDP(7) - MDP(6) - MDP(8)) * t45 + (t44 * MDP(10) - pkin(2) * MDP(7) + MDP(5) - MDP(9)) * qJDD(2); qJDD(2) * MDP(8) - t45 * MDP(9) + (t29 * qJD(2) - t47 + t49) * MDP(10);];
tau  = t1;
