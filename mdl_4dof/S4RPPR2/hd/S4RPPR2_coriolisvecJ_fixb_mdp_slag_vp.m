% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR2
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
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:42
% EndTime: 2019-03-08 18:28:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (123->40), mult. (230->62), div. (0->0), fcn. (110->4), ass. (0->25)
t38 = sin(pkin(6));
t39 = cos(pkin(6));
t60 = qJ(2) * MDP(6) + t38 * MDP(7) + t39 * MDP(8) + MDP(5);
t42 = -pkin(1) - pkin(2);
t58 = t38 * qJ(2);
t55 = qJD(1) * qJ(2);
t54 = qJD(1) * qJD(2);
t52 = t39 * t42 - t58;
t35 = t42 * qJD(1) + qJD(2);
t34 = t39 * t35;
t29 = t34 + (-pkin(3) - t58) * qJD(1);
t31 = t38 * t35 + t39 * t55;
t40 = sin(qJ(4));
t41 = cos(qJ(4));
t51 = t41 * t29 - t40 * t31;
t50 = -t40 * t29 - t41 * t31;
t49 = (-t38 * t55 + t34) * t38 - t31 * t39;
t48 = t38 * t41 + t39 * t40;
t47 = t38 * t40 - t39 * t41;
t45 = t48 * qJD(1);
t44 = t47 * qJD(1);
t37 = -qJD(1) + qJD(4);
t33 = t39 * qJ(2) + t38 * t42;
t32 = -pkin(3) + t52;
t1 = [((t39 * t33 - t38 * t52) * qJD(1) - t49) * qJD(2) * MDP(9) + (((-t32 * t40 - t33 * t41) * t37 - t50) * qJD(4) + (-t48 * t37 + t45) * qJD(2)) * MDP(11) + ((-(t32 * t41 - t33 * t40) * t37 + t51) * qJD(4) + (t47 * t37 - t44) * qJD(2)) * MDP(12) + 0.2e1 * t60 * t54; ((-t48 * qJD(4) + t45) * MDP(11) + (t47 * qJD(4) - t44) * MDP(12)) * t37 + (t49 * MDP(9) - t60 * qJD(1)) * qJD(1); 0; (-t48 * MDP(11) + t47 * MDP(12)) * t54 + (t37 - qJD(4)) * (-t50 * MDP(11) + t51 * MDP(12));];
tauc  = t1;
