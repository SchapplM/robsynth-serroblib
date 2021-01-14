% Calculate Gravitation load on the joints for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [1x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S1P1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S1P1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(1,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [1 1]), ...
  'S1P1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [1x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:56
% EndTime: 2021-01-14 12:21:56
% DurationCPUTime: 0.02s
% Computational Cost: add. (1->1), mult. (1->1), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [-g(3) * MDP(1);];
taug = t1;
