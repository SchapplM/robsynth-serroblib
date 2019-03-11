% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:49
% EndTime: 2019-03-08 18:15:49
% DurationCPUTime: 0.10s
% Computational Cost: add. (47->19), mult. (110->33), div. (0->0), fcn. (68->4), ass. (0->15)
t25 = qJD(3) + qJD(4);
t26 = sin(qJ(4));
t27 = sin(qJ(3));
t42 = t26 * t27;
t45 = t25 * MDP(8) * t42;
t28 = cos(qJ(4));
t43 = t26 * MDP(7) + t28 * MDP(8);
t41 = t27 * t28;
t29 = cos(qJ(3));
t38 = qJD(3) * t29;
t34 = qJD(2) * t45;
t33 = -t26 * t29 - t41;
t32 = t28 * t29 - t42;
t24 = qJD(3) * pkin(3) + t29 * qJD(2);
t1 = [0; (-t27 * MDP(4) - t29 * MDP(5)) * qJD(3) ^ 2 + t25 ^ 2 * (t33 * MDP(7) - t32 * MDP(8)); t34 + ((-t26 * t38 + (-t33 - t41) * t25) * MDP(7) + (t32 * t25 - t28 * t38) * MDP(8)) * qJD(2) + t43 * qJD(4) * (-pkin(3) * t25 - t24); t34 + (-t43 * t38 - t45) * qJD(2) + t43 * t24 * (-qJD(4) + t25);];
tauc  = t1;
