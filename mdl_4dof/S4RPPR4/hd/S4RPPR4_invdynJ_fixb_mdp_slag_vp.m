% Calculate vector of inverse dynamics joint torques for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:55
% EndTime: 2019-12-31 16:38:56
% DurationCPUTime: 0.33s
% Computational Cost: add. (167->69), mult. (266->98), div. (0->0), fcn. (124->8), ass. (0->34)
t94 = 2 * qJD(4);
t59 = qJ(1) + pkin(6);
t57 = sin(t59);
t58 = cos(t59);
t93 = -g(1) * t57 + g(2) * t58 + qJDD(3);
t63 = sin(pkin(6));
t54 = pkin(1) * t63 + qJ(3);
t49 = qJD(3) * qJD(1) + qJDD(1) * t54;
t67 = cos(qJ(4));
t61 = t67 ^ 2;
t65 = sin(qJ(4));
t91 = (t65 ^ 2 - t61) * MDP(9);
t52 = qJD(1) * t54;
t89 = -qJD(1) * t52 + t93;
t69 = qJD(4) ^ 2;
t70 = qJD(1) ^ 2;
t84 = -t69 - t70;
t83 = t67 * MDP(8);
t82 = qJDD(2) - g(3);
t64 = cos(pkin(6));
t56 = -pkin(1) * t64 - pkin(2);
t53 = -pkin(5) + t56;
t80 = qJDD(4) * t53;
t79 = qJDD(4) * t65;
t78 = qJDD(4) * t67;
t76 = -g(1) * t58 - g(2) * t57;
t75 = t56 * MDP(7);
t73 = -t53 * qJDD(1) - t89;
t72 = -t53 * t69 + 0.2e1 * t49 + t76;
t68 = cos(qJ(1));
t66 = sin(qJ(1));
t51 = -t65 * t69 + t78;
t50 = -t67 * t69 - t79;
t1 = [(g(1) * t68 + g(2) * t66) * MDP(3) + t93 * MDP(5) + t76 * MDP(6) + (t49 * t54 + t52 * qJD(3) + qJDD(3) * t56 - g(1) * (-pkin(1) * t66 - pkin(2) * t57 + qJ(3) * t58) - g(2) * (pkin(1) * t68 + pkin(2) * t58 + qJ(3) * t57)) * MDP(7) + t51 * MDP(10) + t50 * MDP(11) + 0.2e1 * (qJD(3) * MDP(6) + qJD(4) * t91) * qJD(1) + (MDP(1) + 0.2e1 * t54 * MDP(6) + t61 * MDP(8) + (t63 ^ 2 + t64 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * MDP(5) + t75) * t56) * qJDD(1) + ((t52 * t94 + t80) * MDP(13) + t72 * MDP(14)) * t67 + (-0.2e1 * qJDD(1) * t67 * MDP(9) + t72 * MDP(13) - MDP(14) * t80 + (-t52 * MDP(14) - qJD(1) * t83) * t94) * t65 + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t66 - g(2) * t68); MDP(13) * t50 - MDP(14) * t51 + (MDP(4) + MDP(7)) * t82; -t70 * MDP(6) + t89 * MDP(7) + (t84 * t65 + t78) * MDP(13) + (t84 * t67 - t79) * MDP(14) + (MDP(5) + t75) * qJDD(1); qJDD(4) * MDP(12) - t70 * t91 + (qJDD(1) * MDP(10) - t73 * MDP(13) - t82 * MDP(14)) * t67 + (-qJDD(1) * MDP(11) - t82 * MDP(13) + t73 * MDP(14) + t70 * t83) * t65;];
tau = t1;
