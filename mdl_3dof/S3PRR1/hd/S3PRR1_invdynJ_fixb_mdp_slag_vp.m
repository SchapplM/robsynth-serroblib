% Calculate vector of inverse dynamics joint torques for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3PRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:01
% EndTime: 2019-03-08 18:04:01
% DurationCPUTime: 0.20s
% Computational Cost: add. (100->51), mult. (158->77), div. (0->0), fcn. (106->6), ass. (0->27)
t48 = qJD(2) + qJD(3);
t67 = t48 ^ 2;
t47 = qJDD(2) + qJDD(3);
t66 = pkin(2) * t47;
t51 = sin(qJ(2));
t65 = qJD(3) * t51;
t53 = cos(qJ(2));
t64 = t53 * qJD(2);
t63 = -qJD(3) + t48;
t62 = qJDD(1) - g(2);
t61 = t51 * qJDD(1);
t60 = qJD(1) * qJD(2);
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t59 = t50 * t53 + t51 * t52;
t58 = -t50 * t51 + t52 * t53;
t39 = qJD(2) * pkin(2) + qJD(1) * t53;
t57 = -qJD(3) * t39 - t61;
t56 = -pkin(2) * qJD(3) * t48 + t57;
t44 = t53 * qJDD(1);
t37 = qJDD(2) * pkin(2) - t51 * t60 + t44;
t49 = qJ(2) + qJ(3);
t45 = sin(t49);
t46 = cos(t49);
t55 = t47 * MDP(5) + (t50 * qJD(1) * t65 + g(1) * t46 + g(2) * t45) * MDP(7) + (g(1) * t45 - g(2) * t46 + t52 * t37) * MDP(6);
t54 = qJD(2) ^ 2;
t1 = [t62 * MDP(1) + (qJDD(2) * t53 - t51 * t54) * MDP(3) + (-qJDD(2) * t51 - t53 * t54) * MDP(4) + (t58 * t47 - t67 * t59) * MDP(6) + (-t59 * t47 - t67 * t58) * MDP(7); qJDD(2) * MDP(2) + (g(1) * t51 - g(2) * t53 + t44) * MDP(3) + (g(1) * t53 - t62 * t51) * MDP(4) + (MDP(6) * t66 + t56 * MDP(7)) * t52 + (t56 * MDP(6) + (-t37 - t66) * MDP(7)) * t50 + ((t59 * t48 - t50 * t64 - t52 * t65) * MDP(6) + (t58 * t48 - t52 * t64) * MDP(7)) * qJD(1) + t55; ((t63 * t39 - t61) * MDP(7) + (t63 * MDP(6) * t51 - MDP(7) * t64) * qJD(1)) * t52 + ((t39 * t48 - t53 * t60 + t57) * MDP(6) + (-qJD(1) * t51 * t48 - t37) * MDP(7)) * t50 + t55;];
tau  = t1;
