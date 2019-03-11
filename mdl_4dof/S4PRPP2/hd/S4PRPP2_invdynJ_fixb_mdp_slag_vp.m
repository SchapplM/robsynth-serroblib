% Calculate vector of inverse dynamics joint torques for
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:01
% EndTime: 2019-03-08 18:19:01
% DurationCPUTime: 0.19s
% Computational Cost: add. (178->78), mult. (319->98), div. (0->0), fcn. (219->6), ass. (0->40)
t84 = sin(qJ(2));
t85 = cos(qJ(2));
t93 = qJD(1) * qJD(2);
t99 = t84 * qJDD(1) + t85 * t93;
t96 = t85 * qJD(1);
t70 = qJD(2) * pkin(2) + t96;
t83 = cos(pkin(5));
t98 = qJD(1) * t84;
t73 = t83 * t98;
t82 = sin(pkin(5));
t57 = t82 * t70 + t73;
t97 = qJDD(2) * pkin(3);
t95 = qJDD(1) - g(2);
t92 = qJD(4) * qJD(2);
t78 = t85 * qJDD(1);
t64 = qJDD(2) * pkin(2) - t84 * t93 + t78;
t52 = t83 * t64 - t99 * t82;
t53 = t82 * t64 + t99 * t83;
t90 = qJDD(2) * qJ(4) + t53;
t89 = qJDD(4) - t52;
t88 = g(1) * t84 - g(2) * t85;
t66 = t82 * t85 + t83 * t84;
t65 = t82 * t84 - t83 * t85;
t56 = t83 * t70 - t82 * t98;
t79 = qJ(2) + pkin(5);
t76 = sin(t79);
t77 = cos(t79);
t87 = g(1) * t76 - g(2) * t77 - t89;
t86 = qJD(2) ^ 2;
t75 = -t83 * pkin(2) - pkin(3);
t74 = t82 * pkin(2) + qJ(4);
t63 = t65 * qJD(1);
t62 = t65 * qJD(2);
t61 = t82 * t96 + t73;
t60 = t66 * qJD(2);
t55 = qJD(2) * qJ(4) + t57;
t54 = -qJD(2) * pkin(3) + qJD(4) - t56;
t51 = t89 - t97;
t50 = t90 + t92;
t1 = [t95 * MDP(1) + (t85 * qJDD(2) - t86 * t84) * MDP(3) + (-qJDD(2) * t84 - t86 * t85) * MDP(4) + (-t52 * t65 + t53 * t66 - t56 * t60 - t57 * t62 - g(2)) * MDP(5) + (-t60 * qJD(2) - t65 * qJDD(2)) * MDP(6) + (-t62 * qJD(2) + t66 * qJDD(2)) * MDP(7) + (t50 * t66 + t51 * t65 + t54 * t60 - t55 * t62 - g(2)) * MDP(8); (t78 + t88) * MDP(3) + (g(1) * t85 - t95 * t84) * MDP(4) + (t56 * t61 + t57 * t63) * MDP(5) + (t61 * qJD(2) + t87) * MDP(6) + (-g(1) * t77 - g(2) * t76 + t63 * qJD(2) + t90 + 0.2e1 * t92) * MDP(7) + (t50 * t74 + t51 * t75 - t54 * t61 - g(1) * (-t76 * pkin(3) + t77 * qJ(4)) - g(2) * (t77 * pkin(3) + t76 * qJ(4)) + (qJD(4) + t63) * t55) * MDP(8) + (MDP(2) + (pkin(3) - t75) * MDP(6) + t74 * MDP(7)) * qJDD(2) + ((t52 * t83 + t53 * t82 + t88) * MDP(5) + t88 * MDP(8)) * pkin(2); (MDP(5) + MDP(8)) * (qJDD(3) - g(3)); -qJDD(2) * MDP(6) - t86 * MDP(7) + (-t55 * qJD(2) - t87 - t97) * MDP(8);];
tau  = t1;
