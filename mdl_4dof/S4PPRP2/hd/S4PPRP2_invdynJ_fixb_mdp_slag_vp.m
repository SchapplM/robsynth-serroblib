% Calculate vector of inverse dynamics joint torques for
% S4PPRP2
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
%   pkin=[a2,a3,a4,d3,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:21
% EndTime: 2019-03-08 18:13:22
% DurationCPUTime: 0.17s
% Computational Cost: add. (137->59), mult. (254->72), div. (0->0), fcn. (184->6), ass. (0->34)
t72 = sin(qJ(3));
t73 = cos(qJ(3));
t84 = qJD(1) * qJD(3);
t92 = qJDD(1) * t72 + t73 * t84;
t71 = cos(pkin(5));
t70 = sin(pkin(5));
t89 = t72 * t70;
t56 = -t73 * t71 + t89;
t52 = t56 * qJD(1);
t91 = qJD(4) + t52;
t90 = MDP(4) + MDP(6);
t88 = t72 * t71;
t87 = qJDD(3) * pkin(3);
t85 = qJDD(1) * t73;
t83 = qJDD(3) * qJ(4);
t82 = t70 * t85 + t92 * t71;
t81 = qJD(1) * t89;
t79 = t92 * t70 - t71 * t85 + t84 * t88;
t78 = -t52 + t81;
t57 = t73 * t70 + t88;
t53 = t57 * qJD(1);
t68 = pkin(5) + qJ(3);
t66 = sin(t68);
t67 = cos(t68);
t76 = g(1) * t66 - g(2) * t67 - t79;
t75 = g(1) * t67 + g(2) * t66 - t82;
t74 = -qJDD(4) + t76;
t55 = t57 * qJD(3);
t54 = t56 * qJD(3);
t51 = qJD(3) * qJ(4) + t53;
t50 = -qJD(3) * pkin(3) + t91;
t48 = qJDD(4) + t79 - t87;
t47 = t83 + (qJD(4) - t81) * qJD(3) + t82;
t1 = [(qJDD(1) - g(2)) * MDP(1) + (-g(2) + (t70 ^ 2 + t71 ^ 2) * qJDD(1)) * MDP(2) + (t47 * t57 + t48 * t56 + t50 * t55 - t51 * t54 - g(2)) * MDP(8) + (-MDP(5) + MDP(7)) * (-t54 * qJD(3) + t57 * qJDD(3)) + t90 * (-t55 * qJD(3) - t56 * qJDD(3)); (MDP(2) + MDP(8)) * (qJDD(2) - g(3)); qJDD(3) * MDP(3) + t76 * MDP(4) + t75 * MDP(5) + (t74 + 0.2e1 * t87) * MDP(6) + (-t75 + 0.2e1 * t83) * MDP(7) + (t47 * qJ(4) - t48 * pkin(3) - t50 * t53 - g(1) * (-t66 * pkin(3) + t67 * qJ(4)) - g(2) * (t67 * pkin(3) + t66 * qJ(4)) + t91 * t51) * MDP(8) + (t78 * MDP(5) + (0.2e1 * qJD(4) - t78) * MDP(7) + t90 * t53) * qJD(3); -qJDD(3) * MDP(6) - qJD(3) ^ 2 * MDP(7) + (-t51 * qJD(3) - t74 - t87) * MDP(8);];
tau  = t1;
