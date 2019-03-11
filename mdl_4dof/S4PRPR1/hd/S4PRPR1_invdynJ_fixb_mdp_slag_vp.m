% Calculate vector of inverse dynamics joint torques for
% S4PRPR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:08
% EndTime: 2019-03-08 18:21:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (154->58), mult. (172->66), div. (0->0), fcn. (86->4), ass. (0->31)
t56 = pkin(6) + qJ(2);
t53 = sin(t56);
t54 = cos(t56);
t87 = g(1) * t53 - g(2) * t54;
t86 = -qJDD(3) + t87;
t80 = pkin(2) * qJDD(2);
t85 = t80 + t86;
t55 = qJDD(2) - qJDD(4);
t61 = -pkin(2) - pkin(3);
t51 = t61 * qJD(2) + qJD(3);
t76 = (qJD(2) * qJD(3));
t68 = qJD(4) * t51 + t76;
t75 = qJD(2) - qJD(4);
t84 = (qJD(4) * t61 + qJD(3)) * t75 + (t55 + qJDD(2)) * qJ(3) + t68;
t47 = t61 * qJDD(2) + qJDD(3);
t78 = qJ(3) * qJD(4);
t83 = -t61 * t55 + (qJD(2) + t75) * t78 - t47;
t81 = t55 * MDP(8);
t77 = qJ(3) * qJDD(2);
t59 = sin(qJ(4));
t60 = cos(qJ(4));
t45 = -t53 * t59 - t54 * t60;
t46 = -t53 * t60 + t54 * t59;
t71 = g(1) * t46 - g(2) * t45;
t70 = -g(1) * t45 - g(2) * t46;
t69 = g(1) * t54 + g(2) * t53;
t66 = -t69 + (2 * t76);
t65 = -t47 + (qJ(3) * t75 + t78) * qJD(2);
t63 = -t51 * t75 - t68 - t77;
t62 = qJD(2) ^ 2;
t1 = [(MDP(1) + MDP(7)) * (qJDD(1) - g(3)); qJDD(2) * MDP(2) + t87 * MDP(3) + t69 * MDP(4) + (0.2e1 * t80 + t86) * MDP(5) + (t66 + 0.2e1 * t77) * MDP(6) + (t85 * pkin(2) + (t66 + t77) * qJ(3)) * MDP(7) + t81 + (t84 * t59 + t83 * t60 - t71) * MDP(9) + (-t83 * t59 + t84 * t60 - t70) * MDP(10); -qJDD(2) * MDP(5) - t62 * MDP(6) + (-t62 * qJ(3) - t85) * MDP(7) + (t59 * MDP(10) - t60 * MDP(9)) * t55 - (t60 * MDP(10) + t59 * MDP(9)) * t75 ^ 2; -t81 + t71 * MDP(9) + t70 * MDP(10) + (t63 * MDP(10) - t65 * MDP(9)) * t60 + (t65 * MDP(10) + t63 * MDP(9)) * t59;];
tau  = t1;
