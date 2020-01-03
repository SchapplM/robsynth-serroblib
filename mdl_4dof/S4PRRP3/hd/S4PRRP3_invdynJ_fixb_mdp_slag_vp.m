% Calculate vector of inverse dynamics joint torques for
% S4PRRP3
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
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:57
% EndTime: 2019-12-31 16:26:58
% DurationCPUTime: 0.37s
% Computational Cost: add. (236->94), mult. (471->131), div. (0->0), fcn. (234->4), ass. (0->48)
t70 = sin(qJ(3));
t67 = t70 ^ 2;
t71 = cos(qJ(3));
t68 = t71 ^ 2;
t99 = (t67 - t68) * MDP(6);
t61 = t71 * pkin(3) + pkin(2);
t98 = t61 * qJDD(2);
t94 = qJ(4) + pkin(5);
t97 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t94 * qJDD(2);
t96 = -0.2e1 * pkin(2);
t95 = g(3) * t71;
t58 = t94 * t70;
t51 = t71 * qJD(1) - qJD(2) * t58;
t90 = qJD(3) * pkin(3);
t50 = t51 + t90;
t93 = t50 - t51;
t91 = MDP(13) * pkin(3);
t89 = t71 * MDP(5);
t88 = pkin(5) * qJDD(3);
t87 = qJDD(2) * pkin(2);
t86 = qJDD(1) - g(3);
t85 = qJDD(2) * t70;
t84 = qJDD(2) * t71;
t82 = qJD(2) * qJD(3);
t83 = t70 * pkin(3) * t82 + qJDD(4);
t59 = t94 * t71;
t81 = t94 * qJD(3);
t66 = pkin(6) + qJ(2);
t62 = sin(t66);
t63 = cos(t66);
t80 = g(1) * t63 + g(2) * t62;
t79 = g(1) * t62 - g(2) * t63;
t78 = qJD(2) * t81;
t53 = t70 * qJD(1) + qJD(2) * t59;
t77 = t50 * t70 - t53 * t71;
t72 = qJD(3) ^ 2;
t75 = -pkin(5) * t72 + t79 + 0.2e1 * t87;
t73 = qJD(2) ^ 2;
t74 = t73 * pkin(2) - qJDD(2) * pkin(5) + t80;
t64 = t71 * qJDD(1);
t57 = qJDD(3) * t71 - t72 * t70;
t56 = qJDD(3) * t70 + t72 * t71;
t55 = -t61 * qJD(2) + qJD(4);
t54 = -t70 * qJD(4) - t71 * t81;
t52 = t71 * qJD(4) - t70 * t81;
t49 = (qJDD(1) - t78) * t70 + t97 * t71;
t48 = qJDD(3) * pkin(3) - t97 * t70 - t71 * t78 + t64;
t1 = [t86 * MDP(1) + t57 * MDP(10) - t56 * MDP(11) + (-t77 * qJD(3) + t48 * t71 + t49 * t70 - g(3)) * MDP(13); qJDD(2) * MDP(2) + t79 * MDP(3) + t67 * qJDD(2) * MDP(5) + t56 * MDP(7) + t57 * MDP(8) + (t49 * t59 + t53 * t52 - t48 * t58 + t50 * t54 - (t83 - t87) * t61 - g(1) * (-t62 * t61 + t63 * t94) - g(2) * (t63 * t61 + t62 * t94)) * MDP(13) - 0.2e1 * t82 * t99 + (t75 * MDP(10) + (t82 * t96 - t88) * MDP(11) + (qJD(2) * t52 - qJD(3) * t50 + qJDD(2) * t59 + t58 * t82 + t49) * MDP(12) + t91 * t98) * t71 + (0.2e1 * MDP(6) * t84 - MDP(10) * t88 - t75 * MDP(11) + (-qJD(2) * t54 + qJDD(2) * t58 - t48) * MDP(12) + (t55 * t91 - t53 * MDP(12) + (MDP(10) * t96 - t59 * MDP(12) + 0.2e1 * t89) * qJD(2)) * qJD(3)) * t70 + (MDP(4) - MDP(12)) * t80; MDP(7) * t85 + MDP(8) * t84 + qJDD(3) * MDP(9) + (t74 * t70 + t64 - t95) * MDP(10) + (-t86 * t70 + t74 * t71) * MDP(11) + (-pkin(3) * t85 + (-t90 + t93) * t71 * qJD(2)) * MDP(12) + (t93 * t53 + (-t95 + t48 + (-qJD(2) * t55 + t80) * t70) * pkin(3)) * MDP(13) + (-t70 * t89 + t99) * t73; (t77 * qJD(2) - t79 + t83 - t98) * MDP(13) + (-t67 - t68) * MDP(12) * t73;];
tau = t1;
