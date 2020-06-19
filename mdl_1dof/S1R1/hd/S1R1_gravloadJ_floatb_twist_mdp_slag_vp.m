% Calculate Gravitation load on the joints for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
% MDP [3x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S1R1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S1R1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [3 1]), ...
  'S1R1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [3x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:57
% EndTime: 2020-06-19 09:12:57
% DurationCPUTime: 0.03s
% Computational Cost: add. (3->3), mult. (6->6), div. (0->0), fcn. (4->2), ass. (0->3)
t4 = cos(qJ(1));
t3 = sin(qJ(1));
t1 = [(g(1) * t3 - g(2) * t4) * MDP(2) + (g(1) * t4 + g(2) * t3) * MDP(3);];
taug = t1;
