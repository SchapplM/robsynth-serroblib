% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (128->49), mult. (340->85), div. (0->0), fcn. (150->2), ass. (0->29)
t64 = -2 * qJD(2) * qJD(3);
t45 = sin(qJ(3));
t43 = t45 ^ 2;
t46 = cos(qJ(3));
t44 = t46 ^ 2;
t63 = -t45 * t46 * MDP(5) + (t43 - t44) * MDP(6);
t48 = qJD(2) ^ 2;
t62 = t48 * pkin(2);
t47 = qJD(3) ^ 2;
t61 = t45 * t47;
t60 = t46 * t47;
t59 = -qJ(4) - pkin(5);
t40 = t59 * t45;
t35 = t46 * qJD(1) + qJD(2) * t40;
t56 = qJD(3) * pkin(3);
t34 = t35 + t56;
t58 = t34 - t35;
t54 = qJD(3) * qJD(1);
t41 = t59 * t46;
t51 = qJD(3) * t59;
t50 = pkin(2) * t64;
t49 = (-t46 * pkin(3) - pkin(2)) * qJD(2);
t38 = -t45 * qJD(4) + t46 * t51;
t36 = t46 * qJD(4) + t45 * t51;
t39 = qJD(4) + t49;
t37 = t45 * qJD(1) - qJD(2) * t41;
t33 = qJD(2) * t38 - t45 * t54;
t32 = qJD(2) * t36 + t46 * t54;
t1 = [(-(t47 * MDP(11)) + (qJD(3) * t37 + t33) * MDP(13)) * t46 + (-(t47 * MDP(10)) + (-qJD(3) * t34 + t32) * MDP(13)) * t45; MDP(7) * t60 - MDP(8) * t61 + (-pkin(5) * t60 + t45 * t50) * MDP(10) + (pkin(5) * t61 + t46 * t50) * MDP(11) + (t32 * t46 - t33 * t45 + (-t34 * t46 - t37 * t45) * qJD(3) + (t36 * t46 - t38 * t45 + (-t40 * t46 + t41 * t45) * qJD(3)) * qJD(2)) * MDP(12) + (-t32 * t41 + t33 * t40 + t34 * t38 + t37 * t36 + (t39 + t49) * t45 * t56) * MDP(13) + t63 * t64; t45 * MDP(10) * t62 + (t58 * t37 + (-qJD(2) * t39 * t45 + t33) * pkin(3)) * MDP(13) + t63 * t48 + ((MDP(11) * t62) + (-t56 + t58) * qJD(2) * MDP(12)) * t46; (-t43 - t44) * MDP(12) * t48 + (-t37 * t46 + (t34 + t56) * t45) * MDP(13) * qJD(2);];
tauc = t1;
