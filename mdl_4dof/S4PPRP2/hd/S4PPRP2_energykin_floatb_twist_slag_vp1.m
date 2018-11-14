% Calculate kinetic energy for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:21
% EndTime: 2018-11-14 13:57:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (302->128), mult. (304->145), div. (0->0), fcn. (138->4), ass. (0->60)
t119 = Icges(5,5) / 0.2e1 - Icges(4,4) / 0.2e1;
t118 = rSges(5,1) + pkin(3);
t117 = rSges(5,3) + qJ(4);
t116 = Icges(5,3) / 0.2e1 + Icges(4,2) / 0.2e1;
t115 = Icges(5,1) / 0.2e1 + Icges(4,1) / 0.2e1;
t110 = Icges(5,4) + Icges(4,5);
t109 = Icges(4,6) - Icges(5,6);
t81 = pkin(5) + qJ(3);
t73 = sin(t81);
t108 = t119 * t73;
t74 = cos(t81);
t107 = t119 * t74;
t83 = cos(pkin(5));
t75 = Icges(3,4) * t83;
t82 = sin(pkin(5));
t106 = Icges(3,1) * t82 / 0.2e1 + t75 / 0.2e1;
t105 = t83 / 0.2e1;
t104 = pkin(2) * t82;
t103 = pkin(2) * t83;
t102 = -pkin(4) - qJ(2);
t101 = -t117 * t74 + t118 * t73;
t100 = t117 * t73 + t118 * t74;
t99 = Icges(3,4) * t82;
t96 = V_base(5) * qJ(2) + V_base(1);
t95 = V_base(4) * qJ(1) + V_base(3);
t94 = qJD(1) + V_base(2);
t93 = t116 * t74 - t108;
t92 = t116 * t73 + t107;
t91 = t115 * t73 - t107;
t90 = t115 * t74 + t108;
t89 = -pkin(1) - t103;
t88 = V_base(6) * pkin(1) + t94;
t87 = qJD(2) + t95;
t86 = V_base(6) * t103 + t88;
t85 = V_base(4) * t104 + t87;
t84 = V_base(5) * pkin(4) + (-qJ(1) - t104) * V_base(6) + t96;
t76 = V_base(6) + qJD(3);
t68 = rSges(3,1) * t83 - t82 * rSges(3,2);
t67 = t82 * rSges(3,1) + rSges(3,2) * t83;
t66 = Icges(3,1) * t83 - t99;
t64 = -Icges(3,2) * t82 + t75;
t63 = Icges(3,2) * t83 + t99;
t60 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t59 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t58 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t57 = V_base(6) * rSges(2,1) + V_base(4) * rSges(2,2) + t94;
t56 = rSges(4,1) * t74 - rSges(4,2) * t73;
t53 = rSges(4,1) * t73 + rSges(4,2) * t74;
t38 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,3) + t95;
t37 = -V_base(5) * rSges(2,2) + V_base(1) + (-qJ(1) - rSges(2,3)) * V_base(6);
t36 = V_base(5) * rSges(3,3) + (-qJ(1) - t67) * V_base(6) + t96;
t35 = V_base(6) * t68 + (-qJ(2) - rSges(3,3)) * V_base(4) + t88;
t34 = V_base(4) * t67 + (-pkin(1) - t68) * V_base(5) + t87;
t33 = V_base(5) * rSges(4,3) - t76 * t53 + t84;
t32 = t76 * t56 + (-rSges(4,3) + t102) * V_base(4) + t86;
t31 = V_base(4) * t53 + (-t56 + t89) * V_base(5) + t85;
t30 = V_base(5) * rSges(5,2) + qJD(4) * t73 - t101 * t76 + t84;
t29 = -qJD(4) * t74 + t100 * t76 + (-rSges(5,2) + t102) * V_base(4) + t86;
t28 = t101 * V_base(4) + (t89 - t100) * V_base(5) + t85;
t1 = m(1) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(2) * (t37 ^ 2 + t38 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t34 ^ 2 + t35 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t31 ^ 2 + t32 ^ 2 + t33 ^ 2) / 0.2e1 + m(5) * (t28 ^ 2 + t29 ^ 2 + t30 ^ 2) / 0.2e1 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + t63 * t105 + t82 * t106 + t93 * t74 + t91 * t73) * V_base(5) ^ 2 + ((Icges(1,1) / 0.2e1 + Icges(2,1) / 0.2e1 - t82 * t64 / 0.2e1 + t66 * t105 + t90 * t74 + t92 * t73) * V_base(4) + (Icges(1,4) + Icges(2,5) + (t106 + t64 / 0.2e1) * t83 + (-t63 / 0.2e1 + t66 / 0.2e1) * t82 + (t91 - t92) * t74 + (t90 - t93) * t73) * V_base(5)) * V_base(4) + ((Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t76 + (t109 * V_base(5) + t110 * V_base(4)) * t74 + (-t109 * V_base(4) + t110 * V_base(5)) * t73) * t76 + ((Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,3) / 0.2e1) * V_base(6) + (Icges(3,5) * t82 + Icges(3,6) * t83 + Icges(1,6) - Icges(2,6)) * V_base(5) + (Icges(3,5) * t83 - Icges(3,6) * t82 - Icges(2,4) + Icges(1,5)) * V_base(4)) * V_base(6);
T  = t1;
