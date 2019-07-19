% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:36
% DurationCPUTime: 0.32s
% Computational Cost: add. (214->73), mult. (327->102), div. (0->0), fcn. (123->4), ass. (0->35)
t86 = MDP(5) + MDP(7);
t67 = -pkin(2) - pkin(3);
t62 = qJD(1) + qJD(2);
t64 = sin(qJ(2));
t84 = pkin(1) * qJD(1);
t74 = t64 * t84;
t54 = qJ(3) * t62 + t74;
t63 = sin(qJ(4));
t85 = t63 * t54;
t66 = cos(qJ(2));
t83 = pkin(1) * qJD(2);
t71 = qJD(1) * t83;
t57 = t66 * t71;
t60 = t62 * qJD(3);
t52 = t57 + t60;
t82 = qJD(2) * t64;
t81 = qJD(4) * t63;
t65 = cos(qJ(4));
t80 = qJD(4) * t65;
t79 = -qJD(4) + t62;
t73 = t66 * t84;
t69 = qJD(3) - t73;
t46 = t62 * t67 + t69;
t70 = t64 * t71;
t78 = t46 * t80 + t52 * t65 + t63 * t70;
t77 = t46 * t81 + t52 * t63 + t54 * t80;
t76 = pkin(1) * t82;
t75 = t66 * t83;
t72 = -pkin(1) * t66 - pkin(2);
t68 = t62 * t73 - t57;
t59 = pkin(1) * t64 + qJ(3);
t58 = -pkin(3) + t72;
t56 = qJD(3) + t75;
t53 = -pkin(2) * t62 + t69;
t1 = [(-t62 * t75 - t57) * MDP(6) + (t56 * t62 + t52) * MDP(8) + (t52 * t59 + t54 * t56 + (qJD(1) * t72 + t53) * t76) * MDP(9) + (-(-qJD(4) * t58 - t56) * t79 * t63 + (-(-qJD(4) * t59 + t76) * t79 - t70) * t65 + t77) * MDP(11) + ((t65 * t56 + t63 * t76) * t79 + ((t58 * t65 - t59 * t63) * t79 - t85) * qJD(4) + t78) * MDP(12) + t86 * (-qJD(1) - t62) * t76; t68 * MDP(6) + (0.2e1 * t60 - t68) * MDP(8) + (t52 * qJ(3) + t54 * qJD(3) + (-t54 * t66 + (-pkin(2) * qJD(2) - t53) * t64) * t84) * MDP(9) + (-(-t63 * qJD(3) + (-qJ(3) * t65 - t63 * t67) * qJD(4)) * t79 + (-t65 * t82 + (-t63 * t66 + t64 * t65) * t79) * t84 + t77) * MDP(11) + (-t85 * qJD(4) + t78 - (-t65 * qJD(3) + (qJ(3) * t63 - t65 * t67) * qJD(4) + (t63 * t64 + t65 * t66) * t84) * t79) * MDP(12) + t86 * (-qJD(2) + t62) * t74; -t62 ^ 2 * MDP(8) + (-t54 * t62 + t70) * MDP(9) - (MDP(11) * t63 + MDP(12) * t65) * t79 ^ 2; (t65 * t70 + (-t46 * t63 - t54 * t65) * t79 - t77) * MDP(11) + (t54 * t81 - (t46 * t65 - t85) * t79 - t78) * MDP(12);];
tauc  = t1;
