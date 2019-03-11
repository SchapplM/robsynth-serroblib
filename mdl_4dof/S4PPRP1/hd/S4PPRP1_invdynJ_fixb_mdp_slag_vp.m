% Calculate vector of inverse dynamics joint torques for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:25
% EndTime: 2019-03-08 18:12:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (88->50), mult. (137->59), div. (0->0), fcn. (94->4), ass. (0->23)
t49 = sin(qJ(3));
t59 = t49 * qJD(2);
t60 = qJDD(3) * pkin(3);
t62 = qJD(3) * t59 - t60;
t61 = sin(pkin(5));
t50 = cos(qJ(3));
t58 = t50 * qJD(2);
t57 = qJDD(3) * qJ(4);
t48 = cos(pkin(5));
t56 = -g(1) * t61 + g(2) * t48;
t39 = qJD(3) * qJ(4) + t59;
t55 = (-qJD(3) * pkin(3) + qJD(4) - t58) * t49 + t39 * t50;
t34 = -t48 * t50 - t61 * t49;
t35 = t48 * t49 - t61 * t50;
t44 = t49 * qJDD(2);
t54 = -g(1) * t34 - g(2) * t35 - t44;
t45 = t50 * qJDD(2);
t53 = g(1) * t35 - g(2) * t34 + t45;
t52 = -qJDD(4) + t53;
t51 = qJD(3) ^ 2;
t33 = qJDD(4) - t45 + t62;
t32 = t57 + t44 + (qJD(4) + t58) * qJD(3);
t1 = [(MDP(1) + MDP(2) + MDP(8)) * (qJDD(1) - g(3)); (qJDD(2) + t56) * MDP(2) + (t55 * qJD(3) + t32 * t49 - t33 * t50 + t56) * MDP(8) + (MDP(4) + MDP(6)) * (t50 * qJDD(3) - t51 * t49) + (-MDP(5) + MDP(7)) * (qJDD(3) * t49 + t51 * t50); qJDD(3) * MDP(3) + t53 * MDP(4) + t54 * MDP(5) + (t52 + 0.2e1 * t60) * MDP(6) + (0.2e1 * qJD(3) * qJD(4) - t54 + 0.2e1 * t57) * MDP(7) + (t32 * qJ(4) + t39 * qJD(4) - t33 * pkin(3) - g(1) * (-pkin(3) * t35 - qJ(4) * t34) - g(2) * (pkin(3) * t34 - qJ(4) * t35) - t55 * qJD(2)) * MDP(8); -qJDD(3) * MDP(6) - t51 * MDP(7) + (-t39 * qJD(3) - t52 + t62) * MDP(8);];
tau  = t1;
