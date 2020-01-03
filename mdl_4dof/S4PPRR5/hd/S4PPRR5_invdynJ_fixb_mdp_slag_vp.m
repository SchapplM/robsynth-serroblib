% Calculate vector of inverse dynamics joint torques for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:54
% EndTime: 2019-12-31 16:19:55
% DurationCPUTime: 0.28s
% Computational Cost: add. (129->56), mult. (283->89), div. (0->0), fcn. (181->6), ass. (0->28)
t73 = qJDD(1) - g(3);
t52 = sin(pkin(6));
t53 = cos(pkin(6));
t65 = g(1) * t52 - g(2) * t53 - qJDD(2);
t54 = sin(qJ(4));
t49 = t54 ^ 2;
t56 = cos(qJ(4));
t84 = (-t56 ^ 2 + t49) * MDP(7);
t82 = 2 * qJD(3);
t70 = qJD(4) * t82;
t71 = qJDD(3) * t56;
t72 = qJDD(3) * t54;
t83 = (t54 * t70 - t71) * MDP(11) + (t56 * t70 + t72) * MDP(12);
t59 = qJD(3) ^ 2;
t79 = t54 * t59;
t77 = qJD(3) * pkin(3);
t76 = t56 * MDP(6);
t69 = -2 * t77;
t68 = g(1) * t53 + g(2) * t52;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t62 = -qJDD(3) * pkin(5) + (qJD(3) * t77) + t65 * t55 - t73 * t57;
t58 = qJD(4) ^ 2;
t61 = -0.2e1 * qJDD(3) * pkin(3) + pkin(5) * t58 + t73 * t55 + t65 * t57;
t60 = (-qJDD(4) * t54 + (-t58 - t59) * t56) * MDP(11) + (-qJDD(4) * t56 + t58 * t54 + t79) * MDP(12);
t46 = -t55 * qJDD(3) - t57 * t59;
t45 = -t57 * qJDD(3) + t55 * t59;
t1 = [t46 * MDP(4) + t45 * MDP(5) + (MDP(1) + MDP(2)) * t73 + t83 * t55 + t60 * t57; -t65 * MDP(2) - t45 * MDP(4) + t46 * MDP(5) + t60 * t55 - t83 * t57; (t49 * MDP(6) + MDP(3)) * qJDD(3) + (-t65 * MDP(4) - t73 * MDP(5)) * t57 + (-t73 * MDP(4) + t65 * MDP(5)) * t55 - 0.2e1 * qJD(4) * t84 * qJD(3) + ((t58 * MDP(8)) - t61 * MDP(11) + (-pkin(5) * MDP(12) + MDP(9)) * qJDD(4) + t69 * MDP(12) * qJD(4)) * t56 + (0.2e1 * MDP(7) * t71 - t58 * MDP(9) + t61 * MDP(12) + (-pkin(5) * MDP(11) + MDP(8)) * qJDD(4) + (t69 * MDP(11) + t76 * t82) * qJD(4)) * t54; -t76 * t79 + t59 * t84 + MDP(8) * t72 + MDP(9) * t71 + qJDD(4) * MDP(10) + (t62 * t54 - t68 * t56) * MDP(11) + (t68 * t54 + t62 * t56) * MDP(12);];
tau = t1;
