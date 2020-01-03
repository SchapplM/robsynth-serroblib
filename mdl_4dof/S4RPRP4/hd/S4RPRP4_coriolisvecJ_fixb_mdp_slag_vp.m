% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:54
% EndTime: 2019-12-31 16:43:55
% DurationCPUTime: 0.27s
% Computational Cost: add. (218->64), mult. (509->102), div. (0->0), fcn. (226->4), ass. (0->36)
t71 = cos(qJ(3));
t62 = sin(pkin(6)) * pkin(1) + pkin(5);
t60 = t62 * qJD(1);
t70 = sin(qJ(3));
t88 = t70 * t60;
t54 = t71 * qJD(2) - t88;
t90 = qJD(4) - t54;
t49 = -qJD(3) * pkin(3) + t90;
t87 = t71 * t60;
t55 = t70 * qJD(2) + t87;
t50 = qJD(3) * qJ(4) + t55;
t66 = t70 ^ 2;
t89 = (-t71 ^ 2 + t66) * MDP(6);
t82 = MDP(10) + MDP(12);
t83 = qJD(2) * qJD(3);
t53 = qJD(3) * t87 + t70 * t83;
t85 = t71 * MDP(5);
t74 = pkin(3) * t70 - qJ(4) * t71;
t56 = t74 * qJD(3) - t70 * qJD(4);
t52 = qJD(1) * t56;
t63 = -cos(pkin(6)) * pkin(1) - pkin(2);
t57 = -t71 * pkin(3) - t70 * qJ(4) + t63;
t51 = qJD(1) * t57;
t61 = qJD(1) * t63;
t81 = MDP(11) - MDP(14);
t80 = t54 + t88;
t79 = 0.2e1 * t51;
t78 = -0.2e1 * t52;
t77 = 0.2e1 * t61;
t75 = t62 * MDP(15) + MDP(13);
t73 = qJD(1) ^ 2;
t72 = qJD(3) ^ 2;
t65 = t71 * t83;
t59 = t74 * qJD(1);
t48 = t65 + (qJD(4) - t88) * qJD(3);
t1 = [(t51 * t56 + t52 * t57) * MDP(15) + (t78 * MDP(12) + (-t82 * t62 + MDP(7)) * t72 + t75 * t48) * t71 + (t78 * MDP(14) + (t81 * t62 - MDP(8)) * t72 + t75 * t53) * t70 + (-0.2e1 * qJD(1) * t89 + (t77 * MDP(11) - t79 * MDP(14) + t75 * t49) * t71 + (t77 * MDP(10) + t79 * MDP(12) + 0.2e1 * qJD(1) * t85 - t75 * t50) * t70) * qJD(3); (t48 * t70 - t53 * t71 + (t49 * t70 + t50 * t71) * qJD(3)) * MDP(15) + (-t82 * t70 - t81 * t71) * t72; (t48 * qJ(4) - t49 * t55 + t90 * t50 - t51 * t59) * MDP(15) - t81 * t65 + (-t70 * t85 + t89) * t73 + (t80 * MDP(11) + (0.2e1 * qJD(4) - t80) * MDP(14) + t82 * t55) * qJD(3) + ((-t61 * MDP(10) - t51 * MDP(12) + t59 * MDP(14)) * t70 + (-t61 * MDP(11) + t59 * MDP(12) + t51 * MDP(14)) * t71) * qJD(1) + (-pkin(3) * MDP(15) - t82) * t53; -t70 * t73 * t71 * MDP(12) + (-t66 * t73 - t72) * MDP(14) + (t51 * t70 * qJD(1) - t50 * qJD(3) + t53) * MDP(15);];
tauc = t1;
