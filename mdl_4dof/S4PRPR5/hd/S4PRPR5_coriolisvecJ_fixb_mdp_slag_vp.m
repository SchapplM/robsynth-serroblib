% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR5
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
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:18
% EndTime: 2019-12-31 16:23:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (110->43), mult. (307->82), div. (0->0), fcn. (204->6), ass. (0->33)
t61 = sin(qJ(4));
t63 = cos(qJ(4));
t78 = (t61 ^ 2 - t63 ^ 2) * MDP(7);
t59 = sin(pkin(7));
t60 = cos(pkin(7));
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t51 = t59 * t64 + t60 * t62;
t65 = qJD(4) ^ 2;
t77 = t51 * t65;
t75 = t63 * MDP(6);
t74 = qJD(1) * t62;
t73 = t61 * qJD(4);
t72 = t63 * MDP(12);
t71 = t63 * qJD(4);
t70 = t64 * qJD(1);
t53 = qJD(2) * pkin(2) + t70;
t42 = t60 * t53 - t59 * t74;
t40 = -qJD(2) * pkin(3) - t42;
t50 = t59 * t62 - t60 * t64;
t48 = t50 * qJD(2);
t45 = qJD(1) * t48;
t69 = -t40 * qJD(2) + t45;
t49 = t50 * qJD(1);
t68 = qJD(2) * (-t60 * pkin(2) - pkin(3)) + t40 - t49;
t46 = t51 * qJD(2);
t44 = qJD(1) * t46;
t54 = t60 * t74;
t47 = t59 * t70 + t54;
t67 = qJD(2) * t47 - (t59 * pkin(2) + pkin(5)) * t65 - t44;
t66 = qJD(2) ^ 2;
t43 = t59 * t53 + t54;
t1 = [(-t42 * t46 - t43 * t48 + t44 * t50 - t45 * t51) * MDP(5) + (t48 * t73 - t63 * t77) * MDP(11) + (t48 * t71 + t61 * t77) * MDP(12) + (-t62 * MDP(3) - t64 * MDP(4)) * t66 + ((-t46 * t63 + t50 * t73) * MDP(11) + (t46 * t61 + t50 * t71) * MDP(12)) * qJD(2); (t42 * t47 + t43 * t49 + (-t44 * t60 - t45 * t59) * pkin(2)) * MDP(5) + (t67 * MDP(11) + t65 * MDP(8)) * t63 + (-t67 * MDP(12) - t65 * MDP(9)) * t61 + (-0.2e1 * qJD(2) * t78 + t68 * t72 + (t68 * MDP(11) + 0.2e1 * qJD(2) * t75) * t61) * qJD(4); (-t61 * MDP(11) - t72) * t65; t66 * t78 + t69 * t72 + (t69 * MDP(11) - t66 * t75) * t61;];
tauc = t1;
