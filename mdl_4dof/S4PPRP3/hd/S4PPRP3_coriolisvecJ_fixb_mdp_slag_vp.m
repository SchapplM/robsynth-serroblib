% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:50
% EndTime: 2019-03-08 18:14:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (36->3), mult. (85->8), div. (0->0), fcn. (46->2), ass. (0->6)
t23 = qJD(3) ^ 2;
t30 = t23 * MDP(5);
t24 = (-MDP(6) * pkin(3) - MDP(4)) * t23;
t22 = cos(qJ(3));
t21 = sin(qJ(3));
t1 = [t30 * t21 + t24 * t22; t24 * t21 - t30 * t22; 0; 0;];
tauc  = t1;
