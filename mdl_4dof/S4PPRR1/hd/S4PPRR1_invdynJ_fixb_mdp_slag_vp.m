% Calculate vector of inverse dynamics joint torques for
% S4PPRR1
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
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:49
% EndTime: 2019-03-08 18:15:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (124->59), mult. (186->88), div. (0->0), fcn. (144->8), ass. (0->32)
t53 = qJD(3) + qJD(4);
t75 = t53 ^ 2;
t55 = qJ(3) + qJ(4);
t74 = cos(t55);
t52 = qJDD(3) + qJDD(4);
t73 = pkin(3) * t52;
t61 = cos(qJ(3));
t72 = qJD(3) * t61;
t59 = sin(qJ(3));
t71 = qJD(4) * t59;
t70 = -qJD(4) + t53;
t69 = t59 * qJDD(2);
t68 = qJD(2) * qJD(3);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t67 = t58 * t61 + t60 * t59;
t66 = -t58 * t59 + t60 * t61;
t47 = qJD(3) * pkin(3) + t61 * qJD(2);
t65 = -qJD(4) * t47 - t69;
t64 = -qJD(4) * pkin(3) * t53 + t65;
t51 = sin(t55);
t56 = sin(pkin(6));
t57 = cos(pkin(6));
t41 = -t56 * t51 - t57 * t74;
t42 = t57 * t51 - t56 * t74;
t50 = t61 * qJDD(2);
t43 = qJDD(3) * pkin(3) - t59 * t68 + t50;
t63 = t52 * MDP(6) + (g(1) * t42 - g(2) * t41 + t60 * t43) * MDP(7) + (t58 * qJD(2) * t71 - g(1) * t41 - g(2) * t42) * MDP(8);
t62 = qJD(3) ^ 2;
t45 = -t56 * t61 + t57 * t59;
t44 = -t56 * t59 - t57 * t61;
t1 = [(MDP(1) + MDP(2)) * (qJDD(1) - g(3)); (-g(1) * t56 + g(2) * t57 + qJDD(2)) * MDP(2) + (t61 * qJDD(3) - t62 * t59) * MDP(4) + (-qJDD(3) * t59 - t62 * t61) * MDP(5) + (t66 * t52 - t67 * t75) * MDP(7) + (-t67 * t52 - t66 * t75) * MDP(8); qJDD(3) * MDP(3) + (g(1) * t45 - g(2) * t44 + t50) * MDP(4) + (-g(1) * t44 - g(2) * t45 - t69) * MDP(5) + (MDP(7) * t73 + t64 * MDP(8)) * t60 + (t64 * MDP(7) + (-t43 - t73) * MDP(8)) * t58 + ((t67 * t53 - t58 * t72 - t60 * t71) * MDP(7) + (t66 * t53 - t60 * t72) * MDP(8)) * qJD(2) + t63; ((t70 * t47 - t69) * MDP(8) + (t70 * MDP(7) * t59 - MDP(8) * t72) * qJD(2)) * t60 + ((t47 * t53 - t61 * t68 + t65) * MDP(7) + (-qJD(2) * t59 * t53 - t43) * MDP(8)) * t58 + t63;];
tau  = t1;
