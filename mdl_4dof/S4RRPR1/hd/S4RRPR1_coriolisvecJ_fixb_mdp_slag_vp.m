% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:46
% EndTime: 2019-01-31 13:16:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (197->59), mult. (494->95), div. (0->0), fcn. (274->6), ass. (0->38)
t63 = sin(qJ(2));
t65 = cos(qJ(2));
t83 = t63 * MDP(5) + t65 * MDP(6);
t59 = qJD(1) + qJD(2);
t58 = qJD(4) + t59;
t82 = qJD(4) - t58;
t60 = sin(pkin(7));
t81 = pkin(2) * t60;
t80 = t60 * t63;
t61 = cos(pkin(7));
t79 = t61 * t63;
t77 = pkin(1) * qJD(1);
t54 = t59 * pkin(2) + t65 * t77;
t72 = t63 * t77;
t42 = t61 * t54 - t60 * t72;
t40 = t59 * pkin(3) + t42;
t64 = cos(qJ(4));
t78 = t64 * t40;
t67 = pkin(1) * (-t60 * t65 - t79);
t49 = qJD(2) * t67;
t45 = qJD(1) * t49;
t66 = pkin(1) * (t61 * t65 - t80);
t51 = qJD(2) * t66;
t46 = qJD(1) * t51;
t62 = sin(qJ(4));
t71 = t64 * t45 - t62 * t46;
t57 = t65 * pkin(1) + pkin(2);
t70 = -pkin(1) * t80 + t61 * t57;
t43 = t60 * t54 + t61 * t72;
t69 = -t62 * t40 - t64 * t43;
t41 = qJD(4) * t62 * t43;
t68 = -t62 * t45 - t64 * t46 + t41;
t56 = t61 * pkin(2) + pkin(3);
t52 = pkin(1) * t79 + t60 * t57;
t50 = qJD(1) * t66;
t48 = qJD(1) * t67;
t47 = pkin(3) + t70;
t1 = [(t42 * t49 + t43 * t51 + t45 * t70 + t46 * t52) * MDP(7) + ((t64 * t49 - t62 * t51) * t58 + t71) * MDP(9) + (-(t62 * t49 + t64 * t51) * t58 + t68) * MDP(10) + (((-t47 * t62 - t52 * t64) * t58 + t69) * MDP(9) + (-(t47 * t64 - t52 * t62) * t58 - t78) * MDP(10)) * qJD(4) + t83 * pkin(1) * qJD(2) * (-qJD(1) - t59); (-t42 * t48 - t43 * t50 + (t45 * t61 + t46 * t60) * pkin(2)) * MDP(7) + (-(t64 * t48 - t62 * t50) * t58 + t71) * MDP(9) + ((t62 * t48 + t64 * t50) * t58 + t68) * MDP(10) + (((-t56 * t62 - t64 * t81) * t58 + t69) * MDP(9) + (-(t56 * t64 - t62 * t81) * t58 - t78) * MDP(10)) * qJD(4) + t83 * t77 * (-qJD(2) + t59); 0; (t82 * t69 + t71) * MDP(9) + (t41 + (-t43 * t58 - t45) * t62 + (-t82 * t40 - t46) * t64) * MDP(10);];
tauc  = t1;
