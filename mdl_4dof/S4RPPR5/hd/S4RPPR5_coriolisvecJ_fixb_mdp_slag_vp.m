% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (108->41), mult. (239->66), div. (0->0), fcn. (98->4), ass. (0->30)
t56 = sin(qJ(4));
t57 = cos(qJ(4));
t83 = (t56 ^ 2 - t57 ^ 2) * MDP(11);
t55 = cos(pkin(6));
t82 = qJ(2) * MDP(6) + t55 * MDP(8) + MDP(5);
t81 = 2 * qJD(1);
t58 = (-pkin(1) - pkin(2));
t48 = qJD(1) * t58 + qJD(2);
t80 = t55 * t48;
t54 = sin(pkin(6));
t79 = t55 * qJ(2) + t54 * t58;
t77 = t54 * qJ(2);
t75 = MDP(16) * t57;
t73 = qJD(2) * t55;
t72 = t57 * MDP(10);
t71 = qJD(1) * qJ(2);
t70 = qJD(2) * t81;
t69 = qJD(4) * t81;
t68 = t54 * t70;
t67 = (-t54 * t71 + t80) * t54 - (t54 * t48 + t55 * t71) * t55;
t66 = t55 * t58 - t77;
t41 = -t80 + (pkin(3) + t77) * qJD(1);
t65 = (t41 - t73) * qJD(1);
t64 = -MDP(15) * t57 + MDP(16) * t56;
t63 = MDP(15) * t56 + t75;
t59 = qJD(4) ^ 2;
t62 = -(-pkin(5) + t79) * t59 + t68;
t61 = qJD(4) * (-qJD(1) * (pkin(3) - t66) - t41 - t73);
t60 = qJD(1) ^ 2;
t1 = [MDP(7) * t68 + ((-t54 * t66 + t55 * t79) * qJD(1) - t67) * qJD(2) * MDP(9) - t69 * t83 + (-(t59 * MDP(12)) + t62 * MDP(15) + t61 * MDP(16)) * t57 + ((t59 * MDP(13)) + t61 * MDP(15) - t62 * MDP(16) + t69 * t72) * t56 + t82 * t70; t64 * t59 * t54 + (0.2e1 * qJD(4) * t55 * t63 + MDP(9) * t67) * qJD(1) + ((-MDP(7) + t64) * t54 - t82) * t60; -t63 * t59; t60 * t83 + t65 * t75 + (MDP(15) * t65 - t60 * t72) * t56;];
tauc = t1;
