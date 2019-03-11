% Calculate vector of inverse dynamics joint torques for
% S3RRP1
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
%   pkin=[a2,a3,d1,d2]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:53
% EndTime: 2019-03-08 18:06:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (172->65), mult. (206->79), div. (0->0), fcn. (82->6), ass. (0->40)
t73 = qJ(1) + qJ(2);
t69 = sin(t73);
t70 = cos(t73);
t96 = -g(1) * t70 - g(2) * t69;
t71 = qJDD(1) + qJDD(2);
t95 = t71 * pkin(2);
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t92 = pkin(1) * qJD(1);
t85 = qJD(2) * t92;
t91 = pkin(1) * qJDD(1);
t94 = t74 * t85 - t76 * t91;
t93 = t74 * t91 + t76 * t85;
t90 = qJD(2) * t74;
t89 = qJD(2) * t76;
t88 = t74 * t92;
t87 = t76 * t92;
t72 = qJD(1) + qJD(2);
t86 = t72 * t90;
t84 = t93 + t96;
t66 = t71 * qJ(3);
t68 = t72 * qJD(3);
t48 = t66 + t68 + t93;
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t83 = g(1) * t75 - g(2) * t77;
t82 = g(1) * t69 - g(2) * t70 - t94;
t81 = -qJDD(3) + t82;
t80 = t72 * t87 - t84;
t79 = t81 + t95;
t78 = -g(1) * (-t69 * pkin(2) + t70 * qJ(3)) - g(2) * (t70 * pkin(2) + t69 * qJ(3));
t67 = t71 * MDP(4);
t63 = -t76 * pkin(1) - pkin(2);
t58 = t74 * pkin(1) + qJ(3);
t53 = pkin(1) * t89 + qJD(3);
t52 = t72 * t88;
t51 = t72 * qJ(3) + t88;
t50 = -t72 * pkin(2) + qJD(3) - t87;
t49 = qJDD(3) + t94 - t95;
t1 = [qJDD(1) * MDP(1) + t83 * MDP(2) + (g(1) * t77 + g(2) * t75) * MDP(3) + t67 + t82 * MDP(5) - t84 * MDP(6) + (-t63 * t71 + t79) * MDP(7) + (t53 * t72 + t58 * t71 + t48 + t96) * MDP(8) + (t48 * t58 + t49 * t63 + t51 * t53 + t78) * MDP(9) + ((t71 * t76 - t86) * MDP(5) + (-t71 * t74 - t72 * t89) * MDP(6) - MDP(7) * t86 + (t50 * t90 + t83) * MDP(9)) * pkin(1); t67 + (t52 + t82) * MDP(5) + t80 * MDP(6) + (t52 + t81 + 0.2e1 * t95) * MDP(7) + (0.2e1 * t66 + 0.2e1 * t68 - t80) * MDP(8) + (t48 * qJ(3) + t51 * qJD(3) - t49 * pkin(2) + (-t50 * t74 - t51 * t76) * t92 + t78) * MDP(9); -t71 * MDP(7) - t72 ^ 2 * MDP(8) + (-t51 * t72 - t79) * MDP(9);];
tau  = t1;
