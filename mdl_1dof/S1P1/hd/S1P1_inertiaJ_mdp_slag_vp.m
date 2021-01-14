% Calculate joint inertia matrix for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [1x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S1P1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S1P1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_inertiaJ_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_inertiaJ_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [1 1]), ...
  'S1P1_inertiaJ_mdp_slag_vp: MDP has to be [1x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:56
% EndTime: 2021-01-14 12:21:56
% DurationCPUTime: 0.02s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [MDP(1);];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
