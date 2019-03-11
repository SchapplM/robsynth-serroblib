% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:01
% DurationCPUTime: 0.20s
% Computational Cost: add. (89->37), mult. (240->62), div. (0->0), fcn. (186->6), ass. (0->23)
t43 = sin(pkin(6));
t44 = cos(pkin(6));
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t58 = -t46 * t43 + t48 * t44;
t42 = qJD(3) + qJD(4);
t57 = qJD(4) - t42;
t40 = t48 * t43 + t46 * t44;
t36 = t40 * qJD(1);
t45 = sin(qJ(4));
t56 = t45 * t36;
t47 = cos(qJ(4));
t54 = t47 * t36;
t35 = t58 * qJD(1);
t32 = qJD(3) * pkin(3) + t35;
t52 = -pkin(3) * t42 - t32;
t37 = t58 * qJD(3);
t33 = qJD(1) * t37;
t38 = t40 * qJD(3);
t34 = qJD(1) * t38;
t51 = -t45 * t33 - t47 * t34;
t31 = qJD(4) * t56;
t1 = [((-t45 * t37 - t47 * t38) * MDP(7) - (t47 * t37 - t45 * t38) * MDP(8) + ((-t40 * t47 - t45 * t58) * MDP(7) - (-t40 * t45 + t47 * t58) * MDP(8)) * qJD(4)) * t42 + (-t38 * MDP(4) - t37 * MDP(5)) * qJD(3); 0; (-(-t45 * t35 - t54) * t42 + t51) * MDP(7) + (t31 - t47 * t33 + t45 * t34 + (t47 * t35 - t56) * t42) * MDP(8) + (t52 * MDP(7) * t45 + (-t36 * MDP(7) + t52 * MDP(8)) * t47) * qJD(4) + (t36 * MDP(4) + t35 * MDP(5) + (-t40 * MDP(4) - MDP(5) * t58) * qJD(1)) * qJD(3); (t51 + t57 * (-t45 * t32 - t54)) * MDP(7) + (t31 + (-t36 * t42 + t34) * t45 + (-t57 * t32 - t33) * t47) * MDP(8);];
tauc  = t1;
