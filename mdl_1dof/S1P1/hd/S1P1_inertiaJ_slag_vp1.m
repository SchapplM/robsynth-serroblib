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
% m [2x1]
%   mass of all robot links (including the base)
% rSges [2x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [2x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S1P1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_inertiaJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_inertiaJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1P1_inertiaJ_slag_vp1: m has to be [2x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [2,3]), ...
  'S1P1_inertiaJ_slag_vp1: rSges has to be [2x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [2 6]), ...
  'S1P1_inertiaJ_slag_vp1: Icges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:51
% EndTime: 2021-01-14 12:21:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [m(2);];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
