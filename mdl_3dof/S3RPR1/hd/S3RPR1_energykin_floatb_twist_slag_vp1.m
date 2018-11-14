% Calculate kinetic energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:21
% EndTime: 2018-11-14 10:14:22
% DurationCPUTime: 0.54s
% Computational Cost: add. (233->111), mult. (348->142), div. (0->0), fcn. (264->4), ass. (0->52)
t99 = Icges(2,4) - Icges(3,5);
t98 = Icges(2,1) + Icges(3,1);
t97 = Icges(2,2) + Icges(3,3);
t74 = sin(qJ(1));
t96 = t99 * t74;
t75 = cos(qJ(1));
t95 = t99 * t75;
t94 = -t97 * t75 - t96;
t93 = t97 * t74 - t95;
t92 = t98 * t74 + t95;
t91 = t98 * t75 - t96;
t90 = Icges(3,4) + Icges(2,5);
t89 = Icges(2,6) - Icges(3,6);
t88 = pkin(2) * t74;
t87 = -pkin(4) - rSges(4,3);
t86 = cos(qJ(3));
t85 = sin(qJ(3));
t43 = -t74 * t85 - t75 * t86;
t83 = Icges(4,4) * t43;
t63 = pkin(1) * t75 + t74 * qJ(2);
t69 = V_base(6) + qJD(1);
t81 = t69 * t63 + V_base(2);
t60 = t74 * pkin(1) - qJ(2) * t75;
t80 = V_base(4) * t60 + V_base(3);
t79 = V_base(5) * pkin(3) + V_base(1);
t76 = qJD(2) * t74 + t79;
t66 = -qJD(3) + t69;
t65 = rSges(2,1) * t75 - t74 * rSges(2,2);
t64 = rSges(3,1) * t75 + t74 * rSges(3,3);
t62 = t74 * rSges(2,1) + rSges(2,2) * t75;
t61 = t74 * rSges(3,1) - rSges(3,3) * t75;
t47 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t46 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t45 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t44 = t74 * t86 - t75 * t85;
t41 = Icges(4,4) * t44;
t39 = V_base(5) * rSges(2,3) - t62 * t69 + t79;
t38 = t65 * t69 + V_base(2) + (-rSges(2,3) - pkin(3)) * V_base(4);
t37 = t62 * V_base(4) - t65 * V_base(5) + V_base(3);
t36 = -rSges(4,1) * t43 + rSges(4,2) * t44;
t35 = rSges(4,1) * t44 + rSges(4,2) * t43;
t34 = -Icges(4,1) * t43 + t41;
t33 = Icges(4,1) * t44 + t83;
t32 = Icges(4,2) * t44 - t83;
t31 = Icges(4,2) * t43 + t41;
t28 = V_base(5) * rSges(3,2) + (-t60 - t61) * t69 + t76;
t27 = -qJD(2) * t75 + t69 * t64 + (-rSges(3,2) - pkin(3)) * V_base(4) + t81;
t26 = t61 * V_base(4) + (-t63 - t64) * V_base(5) + t80;
t25 = -t35 * t66 + t87 * V_base(5) + (-t60 - t88) * t69 + t76;
t24 = t66 * t36 + (pkin(2) * t69 - qJD(2)) * t75 + (-pkin(3) - t87) * V_base(4) + t81;
t23 = (t35 + t88) * V_base(4) + (-pkin(2) * t75 - t36 - t63) * V_base(5) + t80;
t1 = m(1) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t37 ^ 2 + t38 ^ 2 + t39 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((-Icges(4,5) * t44 - Icges(4,6) * t43) * V_base(5) + (Icges(4,5) * t43 - Icges(4,6) * t44) * V_base(4) + Icges(4,3) * t66 / 0.2e1) * t66 + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t69 + (t89 * V_base(5) + t90 * V_base(4)) * t75 + (-t89 * V_base(4) + t90 * V_base(5)) * t74) * t69 + ((t31 * t44 - t33 * t43 + t94 * t74 + t92 * t75) * V_base(5) + (t44 * t32 - t43 * t34 + t93 * t74 + t91 * t75) * V_base(4)) * V_base(4) / 0.2e1 + ((t43 * t31 + t44 * t33 + t92 * t74 - t94 * t75) * V_base(5) + (t32 * t43 + t34 * t44 + t91 * t74 - t93 * t75) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
