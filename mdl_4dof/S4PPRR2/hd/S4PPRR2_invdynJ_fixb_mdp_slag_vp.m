% Calculate vector of inverse dynamics joint torques for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:10
% EndTime: 2019-03-08 18:17:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (174->70), mult. (330->95), div. (0->0), fcn. (270->10), ass. (0->37)
t72 = sin(pkin(6));
t73 = cos(pkin(6));
t75 = sin(qJ(3));
t77 = cos(qJ(3));
t79 = t75 * t72 - t73 * t77;
t71 = qJD(3) + qJD(4);
t93 = qJD(4) - t71;
t69 = qJDD(3) + qJDD(4);
t92 = pkin(3) * t69;
t57 = t72 * t77 + t75 * t73;
t53 = t57 * qJD(1);
t74 = sin(qJ(4));
t91 = t53 * t74;
t76 = cos(qJ(4));
t90 = t53 * t76;
t87 = qJDD(1) * t75;
t86 = qJDD(1) * t77;
t70 = pkin(6) + qJ(3);
t68 = qJ(4) + t70;
t63 = sin(t68);
t64 = cos(t68);
t85 = g(1) * t64 + g(2) * t63 + qJD(4) * t91;
t52 = t79 * qJD(1);
t51 = qJD(3) * pkin(3) - t52;
t84 = -pkin(3) * t71 - t51;
t83 = -t72 * t87 + t73 * t86;
t81 = -t57 * t74 - t76 * t79;
t80 = t57 * t76 - t74 * t79;
t55 = t57 * qJD(3);
t48 = qJDD(3) * pkin(3) - qJD(1) * t55 + t83;
t54 = t79 * qJD(3);
t49 = -qJD(1) * t54 + t57 * qJDD(1);
t78 = g(1) * t63 - g(2) * t64 + t76 * t48 - t74 * t49;
t67 = cos(t70);
t66 = sin(t70);
t65 = t69 * MDP(6);
t1 = [(qJDD(1) - g(2)) * MDP(1) + (-g(2) + (t72 ^ 2 + t73 ^ 2) * qJDD(1)) * MDP(2) + (-qJD(3) * t55 - qJDD(3) * t79) * MDP(4) + (qJD(3) * t54 - qJDD(3) * t57) * MDP(5) + ((-t80 * qJD(4) + t54 * t74 - t55 * t76) * t71 + t81 * t69) * MDP(7) + (-(t81 * qJD(4) - t54 * t76 - t55 * t74) * t71 - t80 * t69) * MDP(8); (qJDD(2) - g(3)) * MDP(2); qJDD(3) * MDP(3) + (g(1) * t66 - g(2) * t67 + t83) * MDP(4) + (g(1) * t67 + g(2) * t66 - t72 * t86 - t73 * t87) * MDP(5) + t65 + (t76 * t92 - (t52 * t74 - t90) * t71 + t78) * MDP(7) + (-t76 * t49 + (-t52 * t76 - t91) * t71 + t85 + (-t92 - t48) * t74) * MDP(8) + (t84 * MDP(7) * t74 + (-t53 * MDP(7) + t84 * MDP(8)) * t76) * qJD(4) + (t53 * MDP(4) - t52 * MDP(5) + (-t57 * MDP(4) + t79 * MDP(5)) * qJD(1)) * qJD(3); t65 + (t78 + t93 * (-t51 * t74 - t90)) * MDP(7) + ((-t53 * t71 - t48) * t74 + (-t93 * t51 - t49) * t76 + t85) * MDP(8);];
tau  = t1;
