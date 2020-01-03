% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:34
% DurationCPUTime: 0.28s
% Computational Cost: add. (161->66), mult. (308->91), div. (0->0), fcn. (107->4), ass. (0->36)
t66 = sin(qJ(4));
t68 = cos(qJ(4));
t102 = t68 * t66 * MDP(10) - (t66 ^ 2 - t68 ^ 2) * MDP(11);
t101 = t66 * MDP(15) + t68 * MDP(16);
t100 = t68 * MDP(15) - t66 * MDP(16);
t98 = -MDP(5) + MDP(7);
t63 = qJD(1) + qJD(2);
t62 = t63 ^ 2;
t97 = -0.2e1 * t102 * qJD(4) * t63;
t69 = cos(qJ(2));
t92 = pkin(1) * qJD(2);
t76 = qJD(1) * t92;
t75 = t69 * t76;
t55 = t63 * qJD(3) + t75;
t67 = sin(qJ(2));
t93 = pkin(1) * qJD(1);
t79 = t67 * t93;
t58 = t63 * qJ(3) + t79;
t89 = qJD(4) * t68;
t96 = t55 * t66 + t58 * t89;
t90 = qJD(4) * t66;
t84 = -MDP(6) + MDP(8);
t82 = qJ(3) * qJD(4);
t81 = t67 * t92;
t80 = t69 * t92;
t77 = -t69 * pkin(1) - pkin(2);
t59 = qJD(3) + t80;
t71 = qJD(4) ^ 2;
t74 = t59 * t63 - (-pkin(6) + t77) * t71;
t61 = t67 * pkin(1) + qJ(3);
t73 = t61 * t63 + t81;
t72 = -t58 * t63 + t67 * t76;
t70 = -pkin(2) - pkin(6);
t57 = -t63 * pkin(2) - t69 * t93 + qJD(3);
t53 = t55 * t68;
t1 = [(t75 + (qJD(3) + t59) * t63) * MDP(8) + (t55 * t61 + t58 * t59 + (t77 * qJD(1) + t57) * t81) * MDP(9) + (t74 * t66 + t73 * t89 + t96) * MDP(15) + (t53 + t74 * t68 + (-t58 - t73) * t90) * MDP(16) + t97 + (-t66 * MDP(12) - t68 * MDP(13)) * t71 + (-MDP(6) * t80 + t98 * t81) * (qJD(1) + t63); (t55 * qJ(3) + t58 * qJD(3)) * MDP(9) + t96 * MDP(15) + (-t58 * t90 + t53) * MDP(16) + ((-MDP(16) * t70 - MDP(13)) * t68 + (-MDP(15) * t70 - MDP(12)) * t66) * t71 + t98 * (qJD(2) - t63) * t79 + (0.2e1 * qJD(3) * MDP(8) + (qJD(3) * t66 + t68 * t82) * MDP(15) + (qJD(3) * t68 - t66 * t82) * MDP(16)) * t63 + (((-pkin(2) * qJD(2) - t57) * MDP(9) - t100 * qJD(4)) * t67 + (-t58 * MDP(9) + t84 * qJD(2) + (-t84 - t101) * t63) * t69) * t93 + t97; -t62 * MDP(8) + t72 * MDP(9) + t101 * (-t62 - t71); t100 * t72 + t102 * t62;];
tauc = t1;
