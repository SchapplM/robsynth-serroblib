% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (60->22), mult. (175->47), div. (0->0), fcn. (86->4), ass. (0->19)
t61 = qJD(3) * pkin(3);
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t60 = t39 * MDP(11) - t37 * MDP(12);
t59 = (t37 ^ 2 - t39 ^ 2) * MDP(7);
t42 = qJD(3) ^ 2;
t53 = MDP(12) * t39;
t46 = MDP(11) * t37 + t53;
t58 = 0.2e1 * qJD(3) * qJD(4) * t46 + t42 * MDP(5);
t55 = t39 * MDP(6);
t41 = qJD(4) ^ 2;
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t45 = t40 * qJD(1) + t38 * qJD(2);
t50 = pkin(5) * t41 + t45 * qJD(3);
t49 = qJD(4) * t61;
t48 = t61 * qJD(3);
t43 = -t42 * MDP(4) - t60 * (t41 + t42);
t1 = [t58 * t38 + t43 * t40; t43 * t38 - t58 * t40; (-t50 * MDP(11) - MDP(12) * t49 + t41 * MDP(8)) * t39 + (-MDP(11) * t49 + t50 * MDP(12) - t41 * MDP(9)) * t37 + (t60 * t45 + (-t46 * pkin(3) + 0.2e1 * t37 * t55 - 0.2e1 * t59) * qJD(4)) * qJD(3); t42 * t59 + t48 * t53 + (t48 * MDP(11) - t42 * t55) * t37;];
tauc = t1;
