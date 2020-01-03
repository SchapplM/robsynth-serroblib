% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->17), mult. (111->27), div. (0->0), fcn. (32->2), ass. (0->10)
t22 = sin(qJ(4));
t23 = cos(qJ(4));
t27 = t23 * MDP(13) - t22 * MDP(14);
t39 = qJ(3) * t27 - t23 * t22 * MDP(8) + (t22 ^ 2 - t23 ^ 2) * MDP(9);
t38 = t22 * MDP(13) + t23 * MDP(14);
t28 = -qJ(3) * MDP(7) - MDP(6);
t26 = qJD(2) ^ 2;
t25 = qJD(4) ^ 2;
t24 = -pkin(2) - pkin(5);
t1 = [-t27 * t25; ((-MDP(14) * t24 - MDP(11)) * t23 + (-MDP(13) * t24 - MDP(10)) * t22) * t25 + (0.2e1 * (-t28 + t38) * qJD(3) + 0.2e1 * t39 * qJD(4)) * qJD(2); t28 * t26 + t38 * (-t25 - t26); -t39 * t26;];
tauc = t1;
