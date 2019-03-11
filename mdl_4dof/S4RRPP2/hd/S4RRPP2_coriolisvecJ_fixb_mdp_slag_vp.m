% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPP2
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
%   see S4RRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:05
% EndTime: 2019-03-08 18:34:05
% DurationCPUTime: 0.15s
% Computational Cost: add. (122->44), mult. (212->52), div. (0->0), fcn. (58->2), ass. (0->27)
t67 = MDP(8) + MDP(11);
t66 = MDP(5) + MDP(7) + MDP(10);
t65 = -pkin(2) - pkin(3);
t48 = qJD(1) + qJD(2);
t49 = sin(qJ(2));
t61 = pkin(1) * qJD(1);
t57 = t49 * t61;
t44 = t48 * qJ(3) + t57;
t50 = cos(qJ(2));
t64 = t44 * t50;
t60 = pkin(1) * qJD(2);
t54 = qJD(1) * t60;
t46 = t50 * t54;
t47 = t48 * qJD(3);
t42 = t46 + t47;
t58 = t50 * t60;
t45 = qJD(3) + t58;
t63 = t42 * (t49 * pkin(1) + qJ(3)) + t44 * t45;
t62 = t42 * qJ(3) + t44 * qJD(3);
t59 = t49 * t60;
t56 = t50 * t61;
t55 = -t50 * pkin(1) - pkin(2);
t52 = qJD(3) - t56;
t51 = t48 * t56 - t46;
t43 = -t48 * pkin(2) + t52;
t37 = t65 * t48 + t52;
t1 = [(-t48 * t58 - t46) * MDP(6) + ((t55 * qJD(1) + t43) * t59 + t63) * MDP(9) + ((t37 + (-pkin(3) + t55) * qJD(1)) * t59 + t63) * MDP(12) + t67 * (t45 * t48 + t42) + t66 * (-qJD(1) - t48) * t59; t51 * MDP(6) + (0.2e1 * t47 - t51) * MDP(8) + ((-t64 + (-pkin(2) * qJD(2) - t43) * t49) * t61 + t62) * MDP(9) + (t52 * t48 + t42) * MDP(11) + ((-t64 + (qJD(2) * t65 - t37) * t49) * t61 + t62) * MDP(12) + t66 * (-qJD(2) + t48) * t57; (MDP(9) + MDP(12)) * (-t44 * t48 + t49 * t54) - t67 * t48 ^ 2; 0;];
tauc  = t1;
