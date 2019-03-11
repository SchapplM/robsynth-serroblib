% Calculate vector of inverse dynamics joint torques for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:11
% DurationCPUTime: 0.22s
% Computational Cost: add. (205->76), mult. (361->104), div. (0->0), fcn. (243->6), ass. (0->39)
t71 = qJ(2) + qJ(3);
t67 = sin(t71);
t68 = cos(t71);
t88 = g(1) * t67 - g(2) * t68;
t70 = qJD(2) + qJD(3);
t73 = sin(qJ(2));
t86 = qJD(1) * t73;
t85 = qJD(3) * t70;
t84 = qJDD(1) - g(2);
t83 = qJD(1) * qJD(2);
t72 = sin(qJ(3));
t81 = qJD(3) * t86;
t60 = t72 * t81;
t82 = g(1) * t68 + g(2) * t67 + t60;
t75 = cos(qJ(2));
t66 = t75 * qJDD(1);
t56 = qJDD(2) * pkin(2) - t73 * t83 + t66;
t74 = cos(qJ(3));
t53 = t74 * t56;
t80 = t53 + t88;
t61 = qJD(2) * pkin(2) + t75 * qJD(1);
t51 = t74 * t61 - t72 * t86;
t58 = t72 * t75 + t74 * t73;
t57 = -t72 * t73 + t74 * t75;
t79 = qJD(3) * t61 + t73 * qJDD(1) + t75 * t83;
t78 = -t79 * MDP(6) - t56 * MDP(7);
t77 = -MDP(6) * t81 - t79 * MDP(7);
t76 = qJD(2) ^ 2;
t69 = qJDD(2) + qJDD(3);
t65 = t69 * MDP(5);
t55 = t57 * qJD(1);
t54 = t58 * qJD(1);
t52 = t72 * t61 + t74 * t86;
t50 = t70 * pkin(3) + t51;
t49 = t70 * t58;
t48 = t70 * t57;
t47 = t72 * t56 + t79 * t74 - t60;
t46 = t69 * pkin(3) - t79 * t72 - t74 * t81 + t53;
t1 = [t84 * MDP(1) + (t75 * qJDD(2) - t76 * t73) * MDP(3) + (-qJDD(2) * t73 - t76 * t75) * MDP(4) + (-t49 * t70 + t57 * t69) * MDP(6) + (-t48 * t70 - t58 * t69) * MDP(7) + (t46 * t57 + t47 * t58 + t52 * t48 - t50 * t49 - g(2)) * MDP(8); qJDD(2) * MDP(2) + (g(1) * t73 - g(2) * t75 + t66) * MDP(3) + (g(1) * t75 - t84 * t73) * MDP(4) + t65 + (t54 * t70 + t80) * MDP(6) + (t55 * t70 + t82) * MDP(7) + (t46 * pkin(3) - t52 * t55 + t50 * t54 - g(1) * (-t73 * pkin(2) - pkin(3) * t67) - g(2) * (t75 * pkin(2) + pkin(3) * t68)) * MDP(8) + ((t69 * MDP(6) - MDP(7) * t85 + (qJD(3) * t52 + t46) * MDP(8)) * pkin(2) + t77) * t74 + ((-MDP(6) * t85 - t69 * MDP(7) + (-qJD(3) * t50 + t47) * MDP(8)) * pkin(2) + t78) * t72; t65 + t80 * MDP(6) + (t51 * t70 + t82) * MDP(7) + t77 * t74 + t78 * t72 + (t70 * MDP(6) + (t50 - t51) * MDP(8)) * t52 + (t46 + t88) * MDP(8) * pkin(3); (qJDD(4) - g(3)) * MDP(8);];
tau  = t1;
