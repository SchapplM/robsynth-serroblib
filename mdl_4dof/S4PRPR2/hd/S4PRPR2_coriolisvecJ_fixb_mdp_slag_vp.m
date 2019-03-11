% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:50
% EndTime: 2019-03-08 18:21:50
% DurationCPUTime: 0.18s
% Computational Cost: add. (121->47), mult. (306->81), div. (0->0), fcn. (228->6), ass. (0->27)
t46 = qJD(2) + qJD(4);
t58 = qJD(4) - t46;
t47 = sin(pkin(6));
t57 = pkin(2) * t47;
t50 = sin(qJ(2));
t56 = qJD(1) * t50;
t48 = cos(pkin(6));
t52 = cos(qJ(2));
t42 = t47 * t52 + t48 * t50;
t37 = t42 * qJD(2);
t35 = qJD(1) * t37;
t41 = -t47 * t50 + t48 * t52;
t39 = t41 * qJD(2);
t36 = qJD(1) * t39;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t55 = -t51 * t35 - t49 * t36;
t44 = qJD(2) * pkin(2) + t52 * qJD(1);
t33 = t48 * t44 - t47 * t56;
t31 = qJD(2) * pkin(3) + t33;
t34 = t47 * t44 + t48 * t56;
t54 = -t49 * t31 - t51 * t34;
t45 = t48 * pkin(2) + pkin(3);
t40 = t41 * qJD(1);
t38 = t42 * qJD(1);
t30 = qJD(4) * t49 * t34;
t1 = [(-t33 * t37 + t34 * t39 - t35 * t41 + t36 * t42) * MDP(5) + (-t50 * MDP(3) - t52 * MDP(4)) * qJD(2) ^ 2 + ((-t51 * t37 - t49 * t39) * MDP(7) - (-t49 * t37 + t51 * t39) * MDP(8) + ((-t41 * t49 - t42 * t51) * MDP(7) - (t41 * t51 - t42 * t49) * MDP(8)) * qJD(4)) * t46; (t33 * t38 - t34 * t40 + (-t35 * t48 + t36 * t47) * pkin(2)) * MDP(5) + (-(-t51 * t38 - t49 * t40) * t46 + t55) * MDP(7) + (t30 - t51 * t36 + t49 * t35 + (-t49 * t38 + t51 * t40) * t46) * MDP(8) + (((-t45 * t49 - t51 * t57) * t46 + t54) * MDP(7) + (-(t45 * t51 - t49 * t57) * t46 - t51 * t31) * MDP(8)) * qJD(4); 0; (t58 * t54 + t55) * MDP(7) + (t30 + (-t34 * t46 + t35) * t49 + (-t58 * t31 - t36) * t51) * MDP(8);];
tauc  = t1;
