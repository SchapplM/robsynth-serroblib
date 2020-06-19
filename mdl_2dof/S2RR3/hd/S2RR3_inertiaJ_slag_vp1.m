% Calculate joint inertia matrix for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_inertiaJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_inertiaJ_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_inertiaJ_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR3_inertiaJ_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:23
% EndTime: 2020-06-19 09:14:24
% DurationCPUTime: 0.14s
% Computational Cost: add. (39->16), mult. (48->22), div. (0->0), fcn. (24->4), ass. (0->12)
t10 = qJ(1) + qJ(2);
t8 = sin(t10);
t9 = cos(t10);
t4 = t9 * rSges(3,1) - t8 * rSges(3,2);
t3 = -t8 * rSges(3,1) - t9 * rSges(3,2);
t12 = cos(qJ(1));
t11 = sin(qJ(1));
t6 = t12 * rSges(2,1) - t11 * rSges(2,2);
t5 = -t11 * rSges(2,1) - t12 * rSges(2,2);
t2 = t12 * pkin(1) + t4;
t1 = -t11 * pkin(1) + t3;
t7 = [Icges(2,3) + Icges(3,3) + m(2) * (t5 ^ 2 + t6 ^ 2) + m(3) * (t1 ^ 2 + t2 ^ 2); m(3) * (t3 * t1 + t4 * t2) + Icges(3,3); m(3) * (t3 ^ 2 + t4 ^ 2) + Icges(3,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
