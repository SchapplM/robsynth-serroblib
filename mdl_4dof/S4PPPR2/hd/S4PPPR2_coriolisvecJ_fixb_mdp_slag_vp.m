% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:18
% EndTime: 2019-03-08 18:10:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (15->4), mult. (42->11), div. (0->0), fcn. (26->4), ass. (0->6)
t13 = qJD(4) ^ 2;
t11 = cos(qJ(4));
t10 = sin(qJ(4));
t9 = cos(pkin(5));
t8 = sin(pkin(5));
t1 = [((t10 * t9 - t11 * t8) * MDP(5) - (-t10 * t8 - t11 * t9) * MDP(6)) * t13; 0; (-t10 * MDP(5) - t11 * MDP(6)) * t13; 0;];
tauc  = t1;
