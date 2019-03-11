% Calculate vector of inverse dynamics joint torques for
% S4RPPR1
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
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:42
% EndTime: 2019-03-08 18:27:42
% DurationCPUTime: 0.27s
% Computational Cost: add. (198->76), mult. (249->97), div. (0->0), fcn. (127->8), ass. (0->33)
t66 = sin(pkin(6));
t56 = pkin(1) * t66 + qJ(3);
t84 = qJDD(1) * t56;
t67 = cos(pkin(6));
t57 = -pkin(1) * t67 - pkin(2);
t86 = t57 * qJDD(1);
t61 = qJDD(1) - qJDD(4);
t83 = t61 * MDP(8);
t82 = qJD(3) * qJD(1);
t81 = qJD(1) - qJD(4);
t55 = -pkin(3) + t57;
t48 = t55 * qJDD(1) + qJDD(3);
t80 = t55 * t61 + t48;
t63 = qJ(1) + pkin(6);
t59 = sin(t63);
t60 = cos(t63);
t68 = sin(qJ(4));
t70 = cos(qJ(4));
t46 = -t59 * t68 - t60 * t70;
t47 = -t59 * t70 + t60 * t68;
t79 = -g(1) * t47 + g(2) * t46;
t78 = g(1) * t46 + g(2) * t47;
t69 = sin(qJ(1));
t71 = cos(qJ(1));
t77 = g(1) * t69 - g(2) * t71;
t49 = t55 * qJD(1) + qJD(3);
t51 = t56 * qJD(1);
t76 = t70 * t49 - t68 * t51;
t75 = -t68 * t49 - t70 * t51;
t74 = g(1) * t59 - g(2) * t60 - qJDD(3);
t50 = t82 + t84;
t73 = qJD(3) * t81 + t56 * t61 + t50;
t1 = [qJDD(1) * MDP(1) + t77 * MDP(2) + (g(1) * t71 + g(2) * t69) * MDP(3) + (t77 + (t66 ^ 2 + t67 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t74 - 0.2e1 * t86) * MDP(5) + (-g(1) * t60 - g(2) * t59 + 0.2e1 * t82 + 0.2e1 * t84) * MDP(6) + (t50 * t56 + t51 * qJD(3) + (qJDD(3) + t86) * t57 - g(1) * (-pkin(1) * t69 - pkin(2) * t59 + qJ(3) * t60) - g(2) * (pkin(1) * t71 + pkin(2) * t60 + qJ(3) * t59)) * MDP(7) + t83 + (-t80 * t70 + t73 * t68 + (-(-t55 * t68 - t56 * t70) * t81 - t75) * qJD(4) + t79) * MDP(9) + (t80 * t68 + t73 * t70 + ((t55 * t70 - t56 * t68) * t81 + t76) * qJD(4) + t78) * MDP(10); (MDP(4) + MDP(7)) * (qJDD(2) - g(3)); -qJD(1) ^ 2 * MDP(6) + (-qJD(1) * t51 - t74) * MDP(7) + (MDP(10) * t68 - MDP(9) * t70) * t61 + (t57 * MDP(7) - MDP(5)) * qJDD(1) - (t70 * MDP(10) + t68 * MDP(9)) * t81 ^ 2; -t83 + (t70 * t48 - t68 * t50 + t75 * t81 - t79) * MDP(9) + (-t68 * t48 - t70 * t50 - t76 * t81 - t78) * MDP(10) + (-t76 * MDP(10) + t75 * MDP(9)) * qJD(4);];
tau  = t1;
