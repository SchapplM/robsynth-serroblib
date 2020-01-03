% Calculate vector of inverse dynamics joint torques for
% S4PRPR4
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
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:02
% EndTime: 2019-12-31 16:22:02
% DurationCPUTime: 0.20s
% Computational Cost: add. (125->53), mult. (195->73), div. (0->0), fcn. (83->4), ass. (0->31)
t44 = pkin(6) + qJ(2);
t42 = sin(t44);
t43 = cos(t44);
t74 = g(1) * t42 - g(2) * t43;
t73 = qJDD(3) - t74;
t66 = pkin(2) * qJDD(2);
t72 = t66 - t73;
t50 = -pkin(2) - pkin(5);
t48 = sin(qJ(4));
t49 = cos(qJ(4));
t70 = t48 * t49;
t46 = t49 ^ 2;
t69 = t48 ^ 2 - t46;
t51 = qJD(4) ^ 2;
t52 = qJD(2) ^ 2;
t68 = -t51 - t52;
t67 = t52 * qJ(3);
t65 = qJDD(1) - g(3);
t64 = qJDD(4) * t48;
t63 = qJDD(4) * t49;
t62 = (qJ(3) * qJDD(2));
t61 = (qJD(2) * qJD(4));
t59 = g(1) * t43 + g(2) * t42;
t57 = (2 * qJ(3) * t61) + qJDD(4) * t50;
t56 = (2 * qJD(2) * qJD(3)) - t59;
t55 = -t50 * qJDD(2) + t67 - t73;
t54 = t56 + (2 * t62);
t53 = -t50 * t51 + t54;
t38 = -t51 * t48 + t63;
t37 = -t51 * t49 - t64;
t1 = [t37 * MDP(13) - t38 * MDP(14) + (MDP(1) + MDP(7)) * t65; qJDD(2) * MDP(2) + t74 * MDP(3) + t59 * MDP(4) + (-0.2e1 * t66 + t73) * MDP(5) + t54 * MDP(6) + (t72 * pkin(2) + (t56 + t62) * qJ(3)) * MDP(7) + (t46 * qJDD(2) - 0.2e1 * t61 * t70) * MDP(8) + 0.2e1 * (-qJDD(2) * t70 + t69 * t61) * MDP(9) + t38 * MDP(10) + t37 * MDP(11) + (t53 * t48 + t57 * t49) * MDP(13) + (-t57 * t48 + t53 * t49) * MDP(14); qJDD(2) * MDP(5) - t52 * MDP(6) + (-t67 - t72) * MDP(7) + (t68 * t48 + t63) * MDP(13) + (t68 * t49 - t64) * MDP(14); qJDD(4) * MDP(12) - t69 * MDP(9) * t52 + (qJDD(2) * MDP(10) - t55 * MDP(13) - t65 * MDP(14)) * t49 + (t49 * t52 * MDP(8) - qJDD(2) * MDP(11) - t65 * MDP(13) + t55 * MDP(14)) * t48;];
tau = t1;
