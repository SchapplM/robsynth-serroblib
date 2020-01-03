% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR7
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
%   see S4PRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (67->37), mult. (183->56), div. (0->0), fcn. (77->4), ass. (0->22)
t38 = sin(qJ(2));
t53 = t38 * qJD(1);
t34 = qJD(2) * qJ(3) + t53;
t46 = -t34 + t53;
t63 = qJD(2) * t46;
t42 = qJD(4) ^ 2;
t43 = qJD(2) ^ 2;
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t60 = t37 * MDP(13) + t39 * MDP(14);
t62 = t60 * (t42 + t43) - MDP(7) * t63;
t44 = t39 * MDP(13) - t37 * MDP(14);
t59 = (t37 ^ 2 - t39 ^ 2) * MDP(9);
t58 = qJD(4) * t46;
t40 = cos(qJ(2));
t50 = t40 * qJD(1);
t48 = t39 * t37 * MDP(8);
t32 = (qJD(3) + t50) * qJD(2);
t47 = -(-pkin(2) - pkin(5)) * t42 + t32;
t45 = qJD(3) - t50;
t33 = -qJD(2) * pkin(2) + t45;
t1 = [(t32 * MDP(7) + (-MDP(3) + MDP(5)) * t43 + (t33 * MDP(7) + 0.2e1 * t44 * qJD(4)) * qJD(2)) * t38 + ((-MDP(4) + MDP(6)) * t43 + t62) * t40; (t32 * qJ(3) + t34 * qJD(3) + (-t33 * t38 - t34 * t40) * qJD(1)) * MDP(7) + (-t42 * MDP(11) - MDP(13) * t58 + t47 * MDP(14)) * t39 + (-t42 * MDP(10) + t47 * MDP(13) + MDP(14) * t58) * t37 + (-pkin(2) * MDP(7) * t53 + 0.2e1 * qJD(3) * MDP(6) + (t44 * qJ(3) - 0.2e1 * t48 + 0.2e1 * t59) * qJD(4) + t60 * t45) * qJD(2); -t43 * MDP(6) - t62; (t48 - t59) * t43 + t44 * t63;];
tauc = t1;
