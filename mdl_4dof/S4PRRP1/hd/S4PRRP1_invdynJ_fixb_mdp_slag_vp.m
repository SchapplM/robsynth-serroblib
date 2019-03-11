% Calculate vector of inverse dynamics joint torques for
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:00
% EndTime: 2019-03-08 18:23:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (207->68), mult. (208->80), div. (0->0), fcn. (82->6), ass. (0->41)
t77 = pkin(6) + qJ(2);
t75 = qJ(3) + t77;
t64 = sin(t75);
t65 = cos(t75);
t100 = g(1) * t65 + g(2) * t64;
t76 = qJDD(2) + qJDD(3);
t99 = t76 * pkin(3);
t80 = sin(qJ(3));
t81 = cos(qJ(3));
t96 = pkin(2) * qJD(2);
t89 = qJD(3) * t96;
t95 = pkin(2) * qJDD(2);
t98 = t80 * t89 - t81 * t95;
t97 = t80 * t95 + t81 * t89;
t94 = qJD(3) * t80;
t93 = qJD(3) * t81;
t92 = t80 * t96;
t91 = t81 * t96;
t78 = qJD(2) + qJD(3);
t90 = t78 * t94;
t88 = -t97 + t100;
t70 = t76 * qJ(4);
t72 = t78 * qJD(4);
t50 = t70 + t72 + t97;
t73 = sin(t77);
t74 = cos(t77);
t87 = g(1) * t73 - g(2) * t74;
t86 = g(1) * t64 - g(2) * t65 - t98;
t85 = -qJDD(4) + t86;
t84 = t78 * t91 + t88;
t83 = t85 + t99;
t82 = -g(1) * (-t64 * pkin(3) + t65 * qJ(4)) - g(2) * (t65 * pkin(3) + t64 * qJ(4));
t71 = t76 * MDP(5);
t67 = -t81 * pkin(2) - pkin(3);
t66 = t80 * pkin(2) + qJ(4);
t57 = pkin(2) * t93 + qJD(4);
t54 = t78 * t92;
t53 = t78 * qJ(4) + t92;
t52 = -t78 * pkin(3) + qJD(4) - t91;
t51 = qJDD(4) + t98 - t99;
t1 = [(MDP(1) + MDP(10)) * (qJDD(1) - g(3)); qJDD(2) * MDP(2) + t87 * MDP(3) + (g(1) * t74 + g(2) * t73) * MDP(4) + t71 + t86 * MDP(6) + t88 * MDP(7) + (-t67 * t76 + t83) * MDP(8) + (t57 * t78 + t66 * t76 - t100 + t50) * MDP(9) + (t50 * t66 + t51 * t67 + t53 * t57 + t82) * MDP(10) + ((t76 * t81 - t90) * MDP(6) + (-t76 * t80 - t78 * t93) * MDP(7) - MDP(8) * t90 + (t52 * t94 + t87) * MDP(10)) * pkin(2); t71 + (t54 + t86) * MDP(6) + t84 * MDP(7) + (t54 + t85 + 0.2e1 * t99) * MDP(8) + (0.2e1 * t70 + 0.2e1 * t72 - t84) * MDP(9) + (t50 * qJ(4) + t53 * qJD(4) - t51 * pkin(3) + (-t52 * t80 - t53 * t81) * t96 + t82) * MDP(10); -t76 * MDP(8) - t78 ^ 2 * MDP(9) + (-t53 * t78 - t83) * MDP(10);];
tau  = t1;
