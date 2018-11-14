% Calculate joint inertia matrix for
% S4PPRP4
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:13
% EndTime: 2018-11-14 14:06:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->14), mult. (50->18), div. (0->0), fcn. (20->2), ass. (0->9)
t10 = -rSges(5,1) - pkin(3);
t9 = m(3) + m(4) + m(5);
t8 = cos(qJ(3));
t7 = sin(qJ(3));
t5 = t8 * rSges(4,1) - t7 * rSges(4,2);
t4 = -t7 * rSges(4,1) - t8 * rSges(4,2);
t2 = -t7 * rSges(5,2) - t10 * t8;
t1 = -t8 * rSges(5,2) + t10 * t7;
t3 = [m(2) + t9; 0; t9; m(4) * t4 + m(5) * t1; m(4) * t5 + m(5) * t2; Icges(4,3) + Icges(5,3) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1) t3(2) t3(4) t3(7); t3(2) t3(3) t3(5) t3(8); t3(4) t3(5) t3(6) t3(9); t3(7) t3(8) t3(9) t3(10);];
Mq  = res;
