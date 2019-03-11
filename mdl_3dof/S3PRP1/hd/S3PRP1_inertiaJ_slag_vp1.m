% Calculate joint inertia matrix for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3PRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_inertiaJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_inertiaJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRP1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:00
% EndTime: 2019-03-08 18:03:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (29->17), mult. (61->23), div. (0->0), fcn. (31->2), ass. (0->9)
t8 = rSges(4,1) + pkin(2);
t7 = rSges(4,3) + qJ(3);
t6 = cos(qJ(2));
t5 = sin(qJ(2));
t4 = t6 * rSges(3,1) - t5 * rSges(3,2);
t3 = -t5 * rSges(3,1) - t6 * rSges(3,2);
t2 = t7 * t5 + t8 * t6;
t1 = -t8 * t5 + t7 * t6;
t9 = [m(2) + m(3) + m(4); m(3) * t4 + m(4) * t2; Icges(3,3) + Icges(4,2) + m(3) * (t3 ^ 2 + t4 ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2); -m(4) * t6; m(4) * (t5 * t1 - t6 * t2); m(4) * (t5 ^ 2 + t6 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t9(1) t9(2) t9(4); t9(2) t9(3) t9(5); t9(4) t9(5) t9(6);];
Mq  = res;
