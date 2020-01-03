% Calculate vector of inverse dynamics joint torques for
% S4PPRR3
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
%   see S4PPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:29
% EndTime: 2019-12-31 16:17:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (111->53), mult. (241->82), div. (0->0), fcn. (157->6), ass. (0->26)
t72 = -2 * qJD(3) * qJD(4);
t59 = cos(pkin(6));
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t82 = sin(pkin(6));
t47 = -t59 * t63 - t82 * t61;
t48 = t59 * t61 - t82 * t63;
t88 = -g(1) * t47 - g(2) * t48 - t61 * qJDD(2);
t87 = g(1) * t48 - g(2) * t47 + t63 * qJDD(2);
t60 = sin(qJ(4));
t56 = t60 ^ 2;
t62 = cos(qJ(4));
t86 = (-t62 ^ 2 + t56) * MDP(7);
t84 = qJD(3) * pkin(3);
t83 = t62 * MDP(6);
t81 = (pkin(5) * qJDD(4));
t78 = qJDD(1) - g(3);
t77 = qJDD(3) * t62;
t64 = qJD(4) ^ 2;
t50 = qJDD(4) * t62 - t64 * t60;
t49 = qJDD(4) * t60 + t64 * t62;
t68 = -2 * t84;
t67 = -qJDD(3) * pkin(5) + (qJD(3) * t84) + t88;
t66 = 0.2e1 * qJDD(3) * pkin(3) - pkin(5) * t64 + t87;
t65 = qJD(3) ^ 2;
t1 = [-t50 * MDP(11) + t49 * MDP(12) + (MDP(1) + MDP(2)) * t78; (-g(1) * t82 + g(2) * t59 + qJDD(2)) * MDP(2) + (qJDD(3) * MDP(4) - (t65 * MDP(5)) + (t60 * t72 + t77) * MDP(11) + (-qJDD(3) * t60 + t62 * t72) * MDP(12)) * t63 + (-t65 * MDP(4) - qJDD(3) * MDP(5) + (-t62 * t65 - t49) * MDP(11) + (t60 * t65 - t50) * MDP(12)) * t61; t87 * MDP(4) + t88 * MDP(5) + t49 * MDP(8) + t50 * MDP(9) + (t56 * MDP(6) + MDP(3)) * qJDD(3) + t86 * t72 + (t66 * MDP(11) + (t68 * qJD(4) - t81) * MDP(12)) * t62 + (0.2e1 * MDP(7) * t77 - MDP(11) * t81 - t66 * MDP(12) + (t68 * MDP(11) + 0.2e1 * qJD(3) * t83) * qJD(4)) * t60; qJDD(4) * MDP(10) + t65 * t86 + (-t78 * MDP(11) + t67 * MDP(12) + qJDD(3) * MDP(9)) * t62 + (t67 * MDP(11) + t78 * MDP(12) + qJDD(3) * MDP(8) - t65 * t83) * t60;];
tau = t1;
