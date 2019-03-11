% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:03
% EndTime: 2019-03-08 18:21:03
% DurationCPUTime: 0.09s
% Computational Cost: add. (49->22), mult. (86->39), div. (0->0), fcn. (26->2), ass. (0->13)
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t28 = t19 * MDP(10);
t32 = t18 * MDP(9) + t28;
t31 = (qJ(3) * MDP(7) + MDP(6));
t20 = -pkin(2) - pkin(3);
t25 = qJD(2) - qJD(4);
t27 = qJD(4) + t25;
t26 = qJD(2) * qJ(3);
t23 = t27 * MDP(9);
t22 = qJD(3) * (qJD(2) + t25);
t16 = t20 * qJD(2) + qJD(3);
t1 = [0; (t18 * t22 + (-(-qJ(3) * t19 - t18 * t20) * t25 + t19 * t26 + t18 * t16) * qJD(4)) * MDP(9) + (t19 * t22 + ((-qJ(3) * t18 + t19 * t20) * t25 - t18 * t26 + t19 * t16) * qJD(4)) * MDP(10) + (2 * t31 * qJD(2) * qJD(3)); -(t31 * qJD(2) ^ 2) - t25 ^ 2 * t32; (-t18 * t23 - t27 * t28) * t16 + (-t32 * qJD(3) + (t27 * MDP(10) * t18 - t19 * t23) * qJ(3)) * qJD(2);];
tauc  = t1;
