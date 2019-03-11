% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:10
% EndTime: 2019-03-08 18:33:10
% DurationCPUTime: 0.17s
% Computational Cost: add. (134->58), mult. (338->81), div. (0->0), fcn. (168->4), ass. (0->34)
t64 = sin(pkin(6));
t67 = cos(qJ(2));
t75 = pkin(1) * qJD(2);
t65 = cos(pkin(6));
t66 = sin(qJ(2));
t79 = t65 * t66;
t52 = (t64 * t67 + t79) * t75;
t48 = qJD(1) * t52;
t80 = t64 * t66;
t78 = t65 * t67;
t63 = qJD(1) + qJD(2);
t76 = pkin(1) * qJD(1);
t72 = t67 * t76;
t55 = t63 * pkin(2) + t72;
t73 = t66 * t76;
t58 = t65 * t73;
t47 = t64 * t55 + t58;
t61 = t67 * pkin(1) + pkin(2);
t77 = pkin(1) * t79 + t64 * t61;
t74 = pkin(1) * t80;
t71 = t64 * t73;
t59 = t75 * t78;
t68 = t65 * t61 - t74;
t56 = qJD(1) * t59;
t49 = -qJD(2) * t71 + t56;
t46 = t65 * t55 - t71;
t62 = t63 * qJD(4);
t53 = (t78 - t80) * t76;
t51 = t64 * t72 + t58;
t50 = -qJD(2) * t74 + qJD(4) + t59;
t45 = t49 + t62;
t44 = t63 * qJ(4) + t47;
t43 = -t63 * pkin(3) + qJD(4) - t46;
t1 = [(-t46 * t52 + t47 * t59 - t48 * t68 + t49 * t77) * MDP(7) - t52 * t63 * MDP(8) + (t50 * t63 + t56 + t62) * MDP(9) + (t45 * (qJ(4) + t77) + t44 * t50 + t48 * (-pkin(3) - t68) + t43 * t52) * MDP(10) + (((-qJD(1) - t63) * MDP(6) - qJD(1) * t64 * MDP(8)) * t67 + (-t47 * t64 * MDP(7) - t63 * MDP(5) + (-t65 * MDP(8) - t64 * MDP(9) - MDP(5)) * qJD(1)) * t66) * t75; (t46 * t51 - t47 * t53 + (-t48 * t65 + t49 * t64) * pkin(2)) * MDP(7) + (t51 * t63 - t48) * MDP(8) + (-t53 * t63 + t49 + 0.2e1 * t62) * MDP(9) + (t45 * (t64 * pkin(2) + qJ(4)) + t48 * (-t65 * pkin(2) - pkin(3)) - t43 * t51 + (qJD(4) - t53) * t44) * MDP(10) + (t66 * MDP(5) + t67 * MDP(6)) * (-qJD(2) + t63) * t76; 0; -t63 ^ 2 * MDP(9) + (-t44 * t63 + t48) * MDP(10);];
tauc  = t1;
