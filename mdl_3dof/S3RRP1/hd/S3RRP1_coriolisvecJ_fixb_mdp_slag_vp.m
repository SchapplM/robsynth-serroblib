% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3RRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:57
% EndTime: 2019-03-08 18:06:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (61->31), mult. (119->44), div. (0->0), fcn. (33->2), ass. (0->17)
t37 = MDP(5) + MDP(7);
t30 = cos(qJ(2));
t35 = pkin(1) * qJD(2);
t31 = qJD(1) * t35;
t26 = t30 * t31;
t28 = qJD(1) + qJD(2);
t27 = t28 * qJD(3);
t22 = t26 + t27;
t36 = pkin(1) * qJD(1);
t29 = sin(qJ(2));
t34 = t29 * t35;
t33 = t30 * t35;
t32 = t29 * t36;
t25 = qJD(3) + t33;
t24 = t28 * qJ(3) + t32;
t23 = -t28 * pkin(2) - t30 * t36 + qJD(3);
t1 = [(-t28 * t33 - t26) * MDP(6) + (t25 * t28 + t22) * MDP(8) + (t22 * (t29 * pkin(1) + qJ(3)) + t24 * t25 + (t23 + (-t30 * pkin(1) - pkin(2)) * qJD(1)) * t34) * MDP(9) + t37 * (-qJD(1) - t28) * t34; -t26 * MDP(6) + (t26 + 0.2e1 * t27) * MDP(8) + (t22 * qJ(3) + t24 * qJD(3)) * MDP(9) + t37 * (-qJD(2) + t28) * t32 + ((-t24 * t30 + (-pkin(2) * qJD(2) - t23) * t29) * MDP(9) + (MDP(6) - MDP(8)) * t30 * t28) * t36; -t28 ^ 2 * MDP(8) + (-t24 * t28 + t29 * t31) * MDP(9);];
tauc  = t1;
