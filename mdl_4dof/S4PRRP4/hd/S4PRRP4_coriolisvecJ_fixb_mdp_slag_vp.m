% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRP4
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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:59
% EndTime: 2019-12-31 16:28:00
% DurationCPUTime: 0.31s
% Computational Cost: add. (159->67), mult. (404->106), div. (0->0), fcn. (167->2), ass. (0->33)
t59 = cos(qJ(3));
t58 = sin(qJ(3));
t75 = qJD(2) * t58;
t69 = pkin(5) * t75;
t50 = t59 * qJD(1) - t69;
t82 = qJD(4) - t50;
t45 = -qJD(3) * pkin(3) + t82;
t68 = t59 * qJD(2) * pkin(5);
t51 = t58 * qJD(1) + t68;
t47 = qJD(3) * qJ(4) + t51;
t56 = t58 ^ 2;
t81 = (-t59 ^ 2 + t56) * MDP(6);
t80 = MDP(10) + MDP(12);
t79 = -2 * pkin(2);
t72 = qJD(1) * qJD(3);
t48 = qJD(3) * t68 + t58 * t72;
t77 = pkin(5) * qJD(3);
t76 = t59 * MDP(5);
t63 = pkin(3) * t58 - qJ(4) * t59;
t44 = t63 * qJD(3) - t58 * qJD(4);
t42 = qJD(2) * t44;
t74 = t58 * MDP(12);
t70 = MDP(11) - MDP(14);
t67 = -0.2e1 * t42;
t55 = t59 * t72;
t66 = t50 * qJD(3) - t55;
t52 = -t59 * pkin(3) - t58 * qJ(4) - pkin(2);
t61 = qJD(2) ^ 2;
t60 = qJD(3) ^ 2;
t49 = t63 * qJD(2);
t46 = qJD(2) * t52;
t43 = t55 + (qJD(4) - t69) * qJD(3);
t1 = [(t43 * t58 - t48 * t59 + (t45 * t58 + t47 * t59) * qJD(3)) * MDP(15) + (-t58 * t80 - t70 * t59) * t60; (t42 * t52 + t46 * t44) * MDP(15) + (t60 * MDP(7) + t67 * MDP(12) + t43 * MDP(13) + (t43 * MDP(15) - t60 * t80) * pkin(5)) * t59 + (-t60 * MDP(8) + t48 * MDP(13) + t67 * MDP(14) + (t48 * MDP(15) + t70 * t60) * pkin(5)) * t58 + ((-t59 * MDP(14) + t74) * t46 + (-0.2e1 * t81 + ((MDP(11) * t79) - t52 * MDP(14)) * t59 + ((MDP(10) * t79) + t52 * MDP(12) + 0.2e1 * t76) * t58) * qJD(2) + (pkin(5) * MDP(15) + MDP(13)) * (t45 * t59 - t47 * t58)) * qJD(3); t66 * MDP(11) + (0.2e1 * qJD(3) * qJD(4) - t66) * MDP(14) + (-t48 * pkin(3) + t43 * qJ(4) - t45 * t51 - t46 * t49 + t82 * t47) * MDP(15) + (-t58 * t76 + t81 + (t58 * MDP(10) + t59 * MDP(11)) * pkin(2)) * t61 + ((t49 * MDP(12) + t46 * MDP(14)) * t59 + (MDP(11) * t77 - t46 * MDP(12) + (t49 - t77) * MDP(14)) * t58) * qJD(2) + t80 * (t51 * qJD(3) - t48); -t61 * t59 * t74 + (-t56 * t61 - t60) * MDP(14) + (-t47 * qJD(3) + t46 * t75 + t48) * MDP(15);];
tauc = t1;
