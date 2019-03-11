% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:00
% EndTime: 2019-03-08 18:19:00
% DurationCPUTime: 0.09s
% Computational Cost: add. (85->39), mult. (225->57), div. (0->0), fcn. (154->4), ass. (0->26)
t46 = cos(pkin(5));
t47 = sin(qJ(2));
t53 = qJD(1) * t47;
t44 = t46 * t53;
t45 = sin(pkin(5));
t48 = cos(qJ(2));
t52 = t48 * qJD(1);
t36 = t45 * t52 + t44;
t40 = t45 * t48 + t46 * t47;
t35 = t40 * qJD(2);
t33 = qJD(1) * t35;
t39 = t45 * t47 - t46 * t48;
t54 = t33 * t39;
t43 = qJD(2) * pkin(2) + t52;
t32 = t45 * t43 + t44;
t51 = t45 * t53;
t31 = t46 * t43 - t51;
t49 = qJD(2) ^ 2;
t42 = t46 * qJD(2) * t52;
t38 = t39 * qJD(1);
t37 = t39 * qJD(2);
t34 = -qJD(2) * t51 + t42;
t30 = t42 + (qJD(4) - t51) * qJD(2);
t29 = qJD(2) * qJ(4) + t32;
t28 = -qJD(2) * pkin(3) + qJD(4) - t31;
t1 = [(-t31 * t35 - t32 * t37 + t34 * t40 + t54) * MDP(5) + (t28 * t35 - t29 * t37 + t30 * t40 + t54) * MDP(8) + (-t47 * MDP(3) - t48 * MDP(4)) * t49 + (-t35 * MDP(6) - t37 * MDP(7)) * qJD(2); (t31 * t36 + t32 * t38 + (-t33 * t46 + t34 * t45) * pkin(2)) * MDP(5) + (t30 * (t45 * pkin(2) + qJ(4)) + t33 * (-t46 * pkin(2) - pkin(3)) - t28 * t36 + (qJD(4) + t38) * t29) * MDP(8) + (t42 + (0.2e1 * qJD(4) + t38 - t51) * qJD(2)) * MDP(7); 0; -t49 * MDP(7) + (t40 * qJD(1) - t29) * MDP(8) * qJD(2);];
tauc  = t1;
