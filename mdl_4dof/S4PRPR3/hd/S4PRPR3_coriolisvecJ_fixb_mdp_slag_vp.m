% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:55
% EndTime: 2019-12-31 16:20:56
% DurationCPUTime: 0.24s
% Computational Cost: add. (122->43), mult. (343->77), div. (0->0), fcn. (224->4), ass. (0->27)
t63 = sin(pkin(7));
t64 = cos(pkin(7));
t72 = (t63 ^ 2 + t64 ^ 2) * qJD(2);
t82 = qJ(3) * t72 * MDP(8);
t81 = MDP(7) * t72;
t65 = sin(qJ(4));
t66 = cos(qJ(4));
t50 = t66 * t63 + t65 * t64;
t47 = t50 * qJD(2);
t79 = pkin(5) + qJ(3);
t77 = qJD(2) * t65;
t76 = qJD(2) * t66;
t73 = t63 * t77;
t58 = -t64 * pkin(3) - pkin(2);
t70 = t65 * t63 - t66 * t64;
t69 = t50 * qJD(3);
t49 = t50 * qJD(4);
t68 = t70 * qJD(3);
t56 = qJD(4) * t64 * t76;
t55 = t79 * t64;
t54 = t79 * t63;
t53 = t58 * qJD(2) + qJD(3);
t48 = t70 * qJD(4);
t46 = t70 * qJD(2);
t43 = qJD(2) * t49;
t42 = -qJD(4) * t73 + t56;
t1 = [(-t49 * MDP(14) + t48 * MDP(15)) * qJD(4); (t42 * t50 - t47 * t48) * MDP(9) + (-t42 * t70 - t50 * t43 + t48 * t46 - t47 * t49) * MDP(10) - t48 * qJD(4) * MDP(11) - t49 * qJD(4) * MDP(12) + (t58 * t43 + t53 * t49 + ((t54 * t65 - t55 * t66) * qJD(4) - t69) * qJD(4)) * MDP(14) + (t58 * t42 - t53 * t48 + ((t54 * t66 + t55 * t65) * qJD(4) + t68) * qJD(4)) * MDP(15) + (0.2e1 * t81 + 0.2e1 * t82) * qJD(3); t56 * MDP(15) + ((t63 * t76 + t64 * t77 + t47) * MDP(14) + (-t46 - t73) * MDP(15)) * qJD(4) + (-t81 - t82) * qJD(2); t47 * t46 * MDP(9) + (-t46 ^ 2 + t47 ^ 2) * MDP(10) + (t56 + (t46 - t73) * qJD(4)) * MDP(11) + (-qJD(2) * t69 - t53 * t47) * MDP(14) + (qJD(2) * t68 + t53 * t46) * MDP(15);];
tauc = t1;
