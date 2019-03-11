% Calculate vector of inverse dynamics joint torques for
% S3RPR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:58
% EndTime: 2019-03-08 18:05:58
% DurationCPUTime: 0.21s
% Computational Cost: add. (121->55), mult. (170->65), div. (0->0), fcn. (86->4), ass. (0->30)
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t83 = g(1) * t54 - g(2) * t56;
t82 = -qJDD(2) + t83;
t76 = pkin(1) * qJDD(1);
t81 = t76 + t82;
t51 = qJDD(1) - qJDD(3);
t57 = -pkin(1) - pkin(2);
t46 = t57 * qJD(1) + qJD(2);
t72 = (qJD(1) * qJD(2));
t64 = qJD(3) * t46 + t72;
t71 = qJD(1) - qJD(3);
t80 = (qJD(3) * t57 + qJD(2)) * t71 + (t51 + qJDD(1)) * qJ(2) + t64;
t45 = t57 * qJDD(1) + qJDD(2);
t74 = qJ(2) * qJD(3);
t79 = -t57 * t51 + (qJD(1) + t71) * t74 - t45;
t77 = t51 * MDP(7);
t73 = qJ(2) * qJDD(1);
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t43 = -t54 * t53 - t56 * t55;
t44 = t53 * t56 - t54 * t55;
t67 = g(1) * t44 - g(2) * t43;
t66 = -g(1) * t43 - g(2) * t44;
t65 = g(1) * t56 + g(2) * t54;
t62 = -t65 + (2 * t72);
t61 = -t45 + (qJ(2) * t71 + t74) * qJD(1);
t59 = -t46 * t71 - t64 - t73;
t58 = qJD(1) ^ 2;
t1 = [qJDD(1) * MDP(1) + t83 * MDP(2) + t65 * MDP(3) + (0.2e1 * t76 + t82) * MDP(4) + (t62 + 0.2e1 * t73) * MDP(5) + (t81 * pkin(1) + (t62 + t73) * qJ(2)) * MDP(6) + t77 + (t80 * t53 + t79 * t55 - t67) * MDP(8) + (-t79 * t53 + t80 * t55 - t66) * MDP(9); -qJDD(1) * MDP(4) - t58 * MDP(5) + (-qJ(2) * t58 - t81) * MDP(6) + (-t55 * MDP(8) + t53 * MDP(9)) * t51 - (t53 * MDP(8) + t55 * MDP(9)) * t71 ^ 2; -t77 + t67 * MDP(8) + t66 * MDP(9) + (-t61 * MDP(8) + t59 * MDP(9)) * t55 + (t59 * MDP(8) + t61 * MDP(9)) * t53;];
tau  = t1;
