% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (88->38), mult. (215->61), div. (0->0), fcn. (91->4), ass. (0->23)
t45 = qJD(5) ^ 2;
t46 = qJD(3) ^ 2;
t41 = sin(qJ(3));
t55 = t41 * qJD(2);
t37 = qJD(3) * qJ(4) + t55;
t51 = -t37 + t55;
t62 = t51 * qJD(3);
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t64 = t40 * MDP(14) + t42 * MDP(15);
t65 = -MDP(8) * t62 + t64 * (t45 + t46);
t63 = qJD(5) * t51;
t61 = (t40 ^ 2 - t42 ^ 2) * MDP(10);
t58 = t42 * MDP(9);
t57 = MDP(14) * t42;
t43 = cos(qJ(3));
t53 = t43 * qJD(2);
t35 = (qJD(4) + t53) * qJD(3);
t52 = -(-pkin(3) - pkin(6)) * t45 + t35;
t50 = qJD(4) - t53;
t47 = -MDP(15) * t40 + t57;
t36 = -qJD(3) * pkin(3) + t50;
t1 = [t47 * t45; (t35 * MDP(8) + (-MDP(4) + MDP(6)) * t46 + (t36 * MDP(8) + 0.2e1 * qJD(5) * t47) * qJD(3)) * t41 + ((-MDP(5) + MDP(7)) * t46 + t65) * t43; (t35 * qJ(4) + t37 * qJD(4) + (-t36 * t41 - t37 * t43) * qJD(2)) * MDP(8) + (-t45 * MDP(12) - MDP(14) * t63 + MDP(15) * t52) * t42 + (-t45 * MDP(11) + MDP(14) * t52 + MDP(15) * t63) * t40 + (-pkin(3) * MDP(8) * t55 + 0.2e1 * qJD(4) * MDP(7) + (qJ(4) * t47 - 0.2e1 * t40 * t58 + 0.2e1 * t61) * qJD(5) + t64 * t50) * qJD(3); -t46 * MDP(7) - t65; -t46 * t61 + t62 * t57 + (-MDP(15) * t62 + t46 * t58) * t40;];
tauc = t1;
