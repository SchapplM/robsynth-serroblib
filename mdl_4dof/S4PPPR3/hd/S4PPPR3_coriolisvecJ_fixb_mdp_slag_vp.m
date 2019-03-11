% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:17
% EndTime: 2019-03-08 18:11:17
% DurationCPUTime: 0.02s
% Computational Cost: add. (21->4), mult. (56->12), div. (0->0), fcn. (40->4), ass. (0->7)
t19 = cos(qJ(4));
t18 = sin(qJ(4));
t17 = cos(pkin(5));
t16 = sin(pkin(5));
t13 = (-t16 * t19 - t17 * t18) * qJD(4);
t12 = (-t16 * t18 + t17 * t19) * qJD(4);
t1 = [(-MDP(5) * t12 - MDP(6) * t13) * qJD(4); (MDP(5) * t13 - MDP(6) * t12) * qJD(4); 0; 0;];
tauc  = t1;
