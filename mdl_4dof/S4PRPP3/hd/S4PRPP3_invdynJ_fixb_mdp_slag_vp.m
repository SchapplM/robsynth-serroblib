% Calculate vector of inverse dynamics joint torques for
% S4PRPP3
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
%   pkin=[a2,a3,a4,d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:04
% EndTime: 2019-03-08 18:20:04
% DurationCPUTime: 0.17s
% Computational Cost: add. (136->67), mult. (195->72), div. (0->0), fcn. (88->2), ass. (0->33)
t87 = MDP(6) + MDP(9);
t66 = -pkin(2) - pkin(3);
t86 = t66 * qJDD(2);
t64 = sin(qJ(2));
t79 = t64 * qJD(1);
t55 = qJD(2) * qJ(3) + t79;
t65 = cos(qJ(2));
t85 = t55 * t65;
t57 = t64 * qJDD(1);
t75 = qJDD(2) * qJ(3);
t78 = t65 * qJD(1);
t47 = t75 + t57 + (qJD(3) + t78) * qJD(2);
t84 = t47 * qJ(3) + t55 * qJD(3);
t83 = t65 * pkin(2) + t64 * qJ(3);
t82 = qJD(2) * t64;
t81 = qJDD(2) * pkin(2);
t80 = t55 * qJD(2);
t77 = MDP(5) + MDP(8);
t74 = t47 * t64 + t65 * t80 - g(2);
t58 = t65 * qJDD(1);
t73 = g(1) * t64 - g(2) * t65 + t58;
t56 = qJD(2) * t79;
t72 = qJDD(3) + t56 - t58;
t71 = qJD(3) - t78;
t70 = -qJDD(3) + t73;
t69 = g(1) * t65 + g(2) * t64 - t57;
t67 = qJD(2) ^ 2;
t60 = t65 * qJ(3);
t54 = -qJD(2) * pkin(2) + t71;
t50 = t66 * qJD(2) + t71;
t48 = t72 - t81;
t46 = t72 + t86;
t1 = [(qJDD(1) - g(2)) * MDP(1) + (-t48 * t65 + t54 * t82 + t74) * MDP(7) + (-t46 * t65 + t50 * t82 + t74) * MDP(10) + (MDP(3) + t77) * (t65 * qJDD(2) - t67 * t64) + (-MDP(4) + t87) * (qJDD(2) * t64 + t67 * t65); qJDD(2) * MDP(2) + t73 * MDP(3) + t69 * MDP(4) + (t70 + 0.2e1 * t81) * MDP(5) + (-t48 * pkin(2) - g(1) * (-t64 * pkin(2) + t60) - g(2) * t83 + (-t54 * t64 - t85) * qJD(1) + t84) * MDP(7) + (t70 - 0.2e1 * t86) * MDP(8) + (t46 * t66 - g(1) * (t66 * t64 + t60) - g(2) * (t65 * pkin(3) + t83) + (-t50 * t64 - t85) * qJD(1) + t84) * MDP(10) + t87 * (0.2e1 * qJD(2) * qJD(3) - t69 + 0.2e1 * t75); -t87 * t67 + (t66 * MDP(10) - pkin(2) * MDP(7) - t77) * qJDD(2) + (MDP(7) + MDP(10)) * (t56 - t70 - t80); (qJDD(4) + g(3)) * MDP(10);];
tau  = t1;
