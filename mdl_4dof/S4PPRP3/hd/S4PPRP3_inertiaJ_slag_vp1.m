% Calculate joint inertia matrix for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:29
% EndTime: 2019-03-08 18:14:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->14), mult. (50->18), div. (0->0), fcn. (20->2), ass. (0->9)
t8 = rSges(5,1) + pkin(3);
t7 = m(3) + m(4) + m(5);
t6 = cos(qJ(3));
t5 = sin(qJ(3));
t4 = t6 * rSges(4,1) - t5 * rSges(4,2);
t3 = -t5 * rSges(4,1) - t6 * rSges(4,2);
t2 = -t5 * rSges(5,2) + t8 * t6;
t1 = -t6 * rSges(5,2) - t8 * t5;
t9 = [m(2) + t7; 0; t7; m(4) * t3 + m(5) * t1; m(4) * t4 + m(5) * t2; Icges(4,3) + Icges(5,3) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1) t9(2) t9(4) t9(7); t9(2) t9(3) t9(5) t9(8); t9(4) t9(5) t9(6) t9(9); t9(7) t9(8) t9(9) t9(10);];
Mq  = res;
