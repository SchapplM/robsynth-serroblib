% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:21
% EndTime: 2019-03-08 18:13:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (61->27), mult. (171->41), div. (0->0), fcn. (116->4), ass. (0->25)
t40 = cos(pkin(5));
t42 = cos(qJ(3));
t39 = sin(pkin(5));
t41 = sin(qJ(3));
t52 = t41 * t39;
t43 = -t42 * t40 + t52;
t32 = t43 * qJD(1);
t53 = qJD(4) + t32;
t50 = MDP(4) + MDP(6);
t51 = t41 * t40;
t48 = qJD(1) * qJD(3);
t46 = t42 * t48;
t31 = t39 * t46 + t48 * t51;
t49 = -MDP(5) + MDP(7);
t47 = qJD(1) * t52;
t45 = -t32 + t47;
t44 = t42 * t39 + t51;
t33 = t44 * qJD(1);
t38 = t40 * t46;
t35 = t44 * qJD(3);
t34 = t43 * qJD(3);
t30 = qJD(3) * qJ(4) + t33;
t29 = -qJD(3) * pkin(3) + t53;
t28 = t38 + (qJD(4) - t47) * qJD(3);
t1 = [(t28 * t44 + t29 * t35 - t30 * t34 + t31 * t43) * MDP(8) + (-t49 * t34 - t50 * t35) * qJD(3); 0; (t28 * qJ(4) - t29 * t33 + t53 * t30) * MDP(8) + t49 * t38 + (t45 * MDP(5) + (0.2e1 * qJD(4) - t45) * MDP(7) + t50 * t33) * qJD(3) + (-pkin(3) * MDP(8) - t50) * t31; -qJD(3) ^ 2 * MDP(7) + (-t30 * qJD(3) + t31) * MDP(8);];
tauc  = t1;
