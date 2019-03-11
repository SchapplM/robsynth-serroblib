% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:11
% DurationCPUTime: 0.14s
% Computational Cost: add. (119->46), mult. (279->77), div. (0->0), fcn. (184->4), ass. (0->30)
t47 = qJD(2) + qJD(3);
t48 = sin(qJ(3));
t49 = sin(qJ(2));
t61 = qJD(1) * qJD(2);
t57 = t49 * t61;
t64 = qJD(1) * t49;
t58 = t48 * t64;
t66 = qJD(3) * t58 + t48 * t57;
t50 = cos(qJ(3));
t65 = t50 * MDP(7);
t51 = cos(qJ(2));
t46 = qJD(2) * pkin(2) + t51 * qJD(1);
t63 = qJD(3) * t46;
t62 = qJD(3) * t50;
t60 = t48 * t63;
t59 = t49 * t62;
t56 = t51 * t61;
t39 = t50 * t46 - t58;
t55 = t48 * t51 + t50 * t49;
t54 = -t48 * t49 + t50 * t51;
t40 = t48 * t46 + t50 * t64;
t53 = t55 * qJD(2);
t42 = t54 * qJD(1);
t41 = t55 * qJD(1);
t38 = t47 * pkin(3) + t39;
t37 = -t55 * qJD(3) - t53;
t36 = t47 * t54;
t35 = -t60 + (-t53 - t59) * qJD(1);
t34 = (t56 + t63) * t50 - t66;
t1 = [(t34 * t55 + t35 * t54 + t40 * t36 + t38 * t37) * MDP(8) + (-t49 * MDP(3) - t51 * MDP(4)) * qJD(2) ^ 2 + (t37 * MDP(6) - t36 * MDP(7)) * t47; (t41 * t47 - t48 * t56 - t50 * t57) * MDP(6) + (t42 * t47 - t50 * t56 + t66) * MDP(7) + (t34 * t48 * pkin(2) + t35 * (t50 * pkin(2) + pkin(3)) - t40 * t42 + t38 * t41) * MDP(8) + (-t40 * MDP(6) - t46 * t65 + ((-t38 * t48 + t40 * t50) * MDP(8) + (-t48 * MDP(6) - t65) * t47) * pkin(2)) * qJD(3); (t40 * t47 - t60) * MDP(6) + (t39 * t47 - t46 * t62 + t66) * MDP(7) + (t35 * pkin(3) + (t38 - t39) * t40) * MDP(8) + (-MDP(6) * t59 + (-t55 * MDP(6) - t51 * t65) * qJD(2)) * qJD(1); 0;];
tauc  = t1;
