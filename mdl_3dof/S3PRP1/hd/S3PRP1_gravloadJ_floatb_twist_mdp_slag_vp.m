% Calculate Gravitation load on the joints for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3PRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:11
% EndTime: 2019-03-08 18:03:11
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->11), mult. (24->15), div. (0->0), fcn. (14->2), ass. (0->4)
t8 = cos(qJ(2));
t7 = sin(qJ(2));
t5 = g(1) * t7 - g(2) * t8;
t1 = [(-MDP(1) - MDP(7)) * g(2); (-g(1) * (-t7 * pkin(2) + t8 * qJ(3)) - g(2) * (t8 * pkin(2) + t7 * qJ(3))) * MDP(7) + (MDP(4) - MDP(6)) * (g(1) * t8 + g(2) * t7) + (MDP(3) + MDP(5)) * t5; -t5 * MDP(7);];
taug  = t1;
