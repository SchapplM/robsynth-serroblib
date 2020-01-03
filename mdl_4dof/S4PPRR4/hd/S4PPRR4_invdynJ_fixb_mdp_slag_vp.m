% Calculate vector of inverse dynamics joint torques for
% S4PPRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:42
% EndTime: 2019-12-31 16:18:43
% DurationCPUTime: 0.30s
% Computational Cost: add. (174->66), mult. (379->103), div. (0->0), fcn. (285->10), ass. (0->39)
t97 = 2 * qJDD(3);
t70 = sin(pkin(6));
t72 = cos(pkin(6));
t96 = g(1) * t72 + g(2) * t70;
t73 = sin(qJ(4));
t67 = t73 ^ 2;
t75 = cos(qJ(4));
t95 = (-t75 ^ 2 + t67) * MDP(7);
t69 = sin(pkin(7));
t71 = cos(pkin(7));
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t60 = t74 * t69 - t76 * t71;
t56 = t60 * qJD(1);
t58 = t60 * qJD(3);
t91 = qJD(3) * pkin(3);
t90 = t75 * MDP(6);
t89 = pkin(5) * qJDD(4);
t54 = -t91 + t56;
t87 = t54 - t56 - t91;
t61 = t76 * t69 + t74 * t71;
t59 = t61 * qJD(3);
t86 = -t59 * qJD(3) - t60 * qJDD(3);
t85 = g(1) * t70 - g(2) * t72 - qJDD(2);
t77 = qJD(4) ^ 2;
t84 = t61 * t77 - t86;
t66 = pkin(7) + qJ(3);
t64 = sin(t66);
t65 = cos(t66);
t83 = g(3) * t64 + t96 * t65;
t82 = -g(3) * t65 + t96 * t64;
t81 = 0.2e1 * qJD(4) * t58 - qJDD(4) * t61;
t80 = -(qJDD(3) * pkin(5)) + qJD(1) * t58 - t54 * qJD(3) - t61 * qJDD(1) + t83;
t57 = t61 * qJD(1);
t79 = pkin(3) * t97 - pkin(5) * t77 - qJD(1) * t59 + t57 * qJD(3) - t60 * qJDD(1) + t82;
t78 = qJD(3) ^ 2;
t63 = qJDD(4) * t75 - t77 * t73;
t62 = qJDD(4) * t73 + t77 * t75;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t69 ^ 2 + t71 ^ 2) * qJDD(1)) * MDP(2) + t86 * MDP(4) + (t58 * qJD(3) - t61 * qJDD(3)) * MDP(5) + (-MDP(11) * t84 + MDP(12) * t81) * t75 + (MDP(11) * t81 + MDP(12) * t84) * t73; t63 * MDP(11) - t62 * MDP(12) - MDP(2) * t85; t82 * MDP(4) + t83 * MDP(5) + t62 * MDP(8) + t63 * MDP(9) + (t67 * MDP(6) + MDP(3)) * qJDD(3) + (-t60 * MDP(4) - t61 * MDP(5)) * qJDD(1) + (t57 * MDP(4) - t56 * MDP(5) - 0.2e1 * qJD(4) * t95 + (-t61 * MDP(4) + t60 * MDP(5)) * qJD(1)) * qJD(3) + (t79 * MDP(11) + (t87 * qJD(4) - t89) * MDP(12)) * t75 + (t75 * MDP(7) * t97 - MDP(11) * t89 - t79 * MDP(12) + (t87 * MDP(11) + 0.2e1 * qJD(3) * t90) * qJD(4)) * t73; qJDD(4) * MDP(10) + t78 * t95 + (-MDP(11) * t85 + MDP(12) * t80 + qJDD(3) * MDP(9)) * t75 + (MDP(11) * t80 + MDP(12) * t85 + qJDD(3) * MDP(8) - t78 * t90) * t73;];
tau = t1;
