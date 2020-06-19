% Calculate joint inertia matrix for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
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
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S1R1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_inertiaJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_inertiaJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1R1_inertiaJ_slag_vp1: m has to be [2x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [2,3]), ...
  'S1R1_inertiaJ_slag_vp1: rSges has to be [2x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [2 6]), ...
  'S1R1_inertiaJ_slag_vp1: Icges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:53
% EndTime: 2020-06-19 09:12:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (4->4), mult. (10->7), div. (0->0), fcn. (4->2), ass. (0->5)
t4 = cos(qJ(1));
t3 = sin(qJ(1));
t2 = rSges(2,1) * t4 - rSges(2,2) * t3;
t1 = -rSges(2,1) * t3 - rSges(2,2) * t4;
t5 = [m(2) * (t1 ^ 2 + t2 ^ 2) + Icges(2,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t5(1);];
Mq = res;
