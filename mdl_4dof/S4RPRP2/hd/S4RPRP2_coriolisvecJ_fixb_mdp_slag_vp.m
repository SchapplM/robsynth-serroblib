% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:53
% EndTime: 2019-03-08 18:30:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (131->39), mult. (219->58), div. (0->0), fcn. (77->2), ass. (0->25)
t45 = -pkin(1) - pkin(2);
t40 = t45 * qJD(1) + qJD(2);
t43 = sin(qJ(3));
t58 = t43 * t40;
t50 = qJD(3) * t58;
t54 = t43 * qJD(2);
t44 = cos(qJ(3));
t55 = qJD(3) * t44;
t32 = -t50 + (-qJ(2) * t55 - t54) * qJD(1);
t52 = qJ(2) * qJD(1);
t37 = t44 * t52 + t58;
t51 = qJD(1) - qJD(3);
t59 = -t37 * t51 + t32;
t53 = t44 * qJD(2);
t57 = qJD(1) * t53 + t40 * t55;
t56 = t43 * MDP(8);
t49 = t43 * t52;
t36 = t44 * t40 - t49;
t48 = t44 * qJ(2) + t43 * t45;
t47 = -t43 * qJ(2) + t44 * t45;
t31 = -qJD(3) * t49 + t57;
t35 = -t48 * qJD(3) - t54;
t34 = t47 * qJD(3) + t53;
t33 = -pkin(3) * t51 + t36;
t1 = [(-t35 * t51 + t50) * MDP(8) + (t34 * t51 + t57) * MDP(9) + (t31 * t48 + t37 * t34 + t32 * (-pkin(3) + t47) + t33 * t35) * MDP(10) + (((2 * MDP(5)) + t56) * qJD(2) + (0.2e1 * qJD(2) * MDP(6) + (MDP(8) * t44 - MDP(9) * t43) * qJD(3)) * qJ(2)) * qJD(1); (-qJ(2) * MDP(6) - MDP(5)) * qJD(1) ^ 2 - (t44 * MDP(9) + t56) * t51 ^ 2 + (t59 * t44 + (t51 * t33 + t31) * t43) * MDP(10); t59 * MDP(8) + (-t36 * t51 - t31) * MDP(9) + (t32 * pkin(3) + (t33 - t36) * t37) * MDP(10); 0;];
tauc  = t1;
