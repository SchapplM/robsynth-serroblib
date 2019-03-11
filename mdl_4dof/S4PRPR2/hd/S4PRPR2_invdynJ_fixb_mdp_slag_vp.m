% Calculate vector of inverse dynamics joint torques for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:55
% EndTime: 2019-03-08 18:21:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (222->78), mult. (409->112), div. (0->0), fcn. (310->8), ass. (0->44)
t78 = qJD(2) + qJD(4);
t99 = qJD(4) - t78;
t79 = sin(pkin(6));
t98 = pkin(2) * t79;
t82 = sin(qJ(2));
t97 = qJD(1) * t82;
t96 = qJDD(1) - g(2);
t95 = qJD(1) * qJD(2);
t84 = cos(qJ(2));
t67 = qJD(2) * pkin(2) + qJD(1) * t84;
t80 = cos(pkin(6));
t57 = t67 * t79 + t80 * t97;
t76 = qJ(2) + pkin(6) + qJ(4);
t71 = sin(t76);
t72 = cos(t76);
t81 = sin(qJ(4));
t94 = qJD(4) * t81 * t57 + g(1) * t72 + g(2) * t71;
t93 = g(1) * t82 - g(2) * t84;
t56 = t80 * t67 - t79 * t97;
t55 = qJD(2) * pkin(3) + t56;
t83 = cos(qJ(4));
t92 = -t81 * t55 - t83 * t57;
t64 = -t79 * t82 + t80 * t84;
t65 = t79 * t84 + t80 * t82;
t91 = t64 * t83 - t65 * t81;
t90 = t64 * t81 + t65 * t83;
t73 = pkin(2) * t80 + pkin(3);
t89 = t73 * t81 + t83 * t98;
t88 = t73 * t83 - t81 * t98;
t75 = t84 * qJDD(1);
t63 = qJDD(2) * pkin(2) - t82 * t95 + t75;
t86 = qJDD(1) * t82 + t84 * t95;
t52 = t80 * t63 - t86 * t79;
t51 = qJDD(2) * pkin(3) + t52;
t53 = t63 * t79 + t86 * t80;
t87 = g(1) * t71 - g(2) * t72 + t83 * t51 - t81 * t53;
t85 = qJD(2) ^ 2;
t77 = qJDD(2) + qJDD(4);
t74 = t77 * MDP(6);
t62 = t64 * qJD(1);
t61 = t64 * qJD(2);
t60 = t65 * qJD(1);
t59 = t65 * qJD(2);
t1 = [t96 * MDP(1) + (qJDD(2) * t84 - t82 * t85) * MDP(3) + (-qJDD(2) * t82 - t84 * t85) * MDP(4) + (t52 * t64 + t53 * t65 - t56 * t59 + t57 * t61 - g(2)) * MDP(5) + ((-t90 * qJD(4) - t59 * t83 - t61 * t81) * t78 + t91 * t77) * MDP(7) + (-(t91 * qJD(4) - t59 * t81 + t61 * t83) * t78 - t90 * t77) * MDP(8); qJDD(2) * MDP(2) + (t75 + t93) * MDP(3) + (g(1) * t84 - t96 * t82) * MDP(4) + (t56 * t60 - t57 * t62 + (t52 * t80 + t53 * t79 + t93) * pkin(2)) * MDP(5) + t74 + (t88 * t77 - (-t60 * t83 - t62 * t81) * t78 + t87) * MDP(7) + (-t89 * t77 - t83 * t53 - t81 * t51 + (-t60 * t81 + t62 * t83) * t78 + t94) * MDP(8) + ((-t89 * t78 + t92) * MDP(7) + (-t83 * t55 - t88 * t78) * MDP(8)) * qJD(4); (qJDD(3) - g(3)) * MDP(5); t74 + (t99 * t92 + t87) * MDP(7) + ((-t57 * t78 - t51) * t81 + (-t99 * t55 - t53) * t83 + t94) * MDP(8);];
tau  = t1;
