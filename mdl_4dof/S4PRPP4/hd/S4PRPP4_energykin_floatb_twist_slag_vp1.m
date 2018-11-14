% Calculate kinetic energy for
% S4PRPP4
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
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:17
% EndTime: 2018-11-14 14:09:18
% DurationCPUTime: 0.43s
% Computational Cost: add. (311->128), mult. (304->147), div. (0->0), fcn. (138->4), ass. (0->59)
t118 = Icges(5,5) / 0.2e1 - Icges(4,4) / 0.2e1;
t117 = rSges(5,1) + pkin(3);
t116 = rSges(5,3) + qJ(4);
t115 = Icges(5,3) / 0.2e1 + Icges(4,2) / 0.2e1;
t114 = Icges(5,1) / 0.2e1 + Icges(4,1) / 0.2e1;
t109 = -Icges(5,4) - Icges(4,5);
t108 = Icges(4,6) - Icges(5,6);
t81 = qJ(2) + pkin(5);
t73 = sin(t81);
t107 = t118 * t73;
t74 = cos(t81);
t106 = t118 * t74;
t83 = cos(qJ(2));
t78 = Icges(3,4) * t83;
t82 = sin(qJ(2));
t105 = Icges(3,1) * t82 / 0.2e1 + t78 / 0.2e1;
t104 = t83 / 0.2e1;
t103 = pkin(2) * t82;
t102 = pkin(2) * t83;
t101 = -pkin(4) - qJ(3);
t100 = -t116 * t74 + t117 * t73;
t99 = t116 * t73 + t117 * t74;
t98 = Icges(3,4) * t82;
t95 = V_base(6) * qJ(1) + V_base(2);
t94 = V_base(4) * pkin(1) + V_base(3);
t93 = qJD(1) + V_base(1);
t92 = t115 * t74 - t107;
t91 = t115 * t73 + t106;
t90 = t114 * t73 - t106;
t89 = t114 * t74 + t107;
t88 = -qJ(1) - t103;
t87 = V_base(4) * pkin(4) + t95;
t86 = V_base(4) * t102 - qJD(3) + t94;
t75 = V_base(6) - qJD(2);
t85 = V_base(4) * qJ(3) + t75 * t103 + t87;
t84 = -V_base(6) * pkin(1) + t93;
t68 = rSges(3,1) * t83 - t82 * rSges(3,2);
t67 = t82 * rSges(3,1) + rSges(3,2) * t83;
t66 = Icges(3,1) * t83 - t98;
t64 = -Icges(3,2) * t82 + t78;
t63 = Icges(3,2) * t83 + t98;
t60 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t59 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t58 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t57 = -V_base(6) * rSges(2,1) + V_base(5) * rSges(2,2) + t93;
t56 = rSges(4,1) * t74 - rSges(4,2) * t73;
t53 = rSges(4,1) * t73 + rSges(4,2) * t74;
t38 = V_base(4) * rSges(2,1) + V_base(3) + (-rSges(2,3) - qJ(1)) * V_base(5);
t37 = -V_base(4) * rSges(2,2) + V_base(6) * rSges(2,3) + t95;
t36 = V_base(4) * rSges(3,3) + t67 * t75 + t87;
t35 = -t68 * t75 + (-pkin(4) - rSges(3,3)) * V_base(5) + t84;
t34 = t68 * V_base(4) + (-qJ(1) - t67) * V_base(5) + t94;
t33 = V_base(4) * rSges(4,3) + t53 * t75 + t85;
t32 = (-t56 - t102) * t75 + (-rSges(4,3) + t101) * V_base(5) + t84;
t31 = t56 * V_base(4) + (-t53 + t88) * V_base(5) + t86;
t30 = V_base(4) * rSges(5,2) + qJD(4) * t73 + t100 * t75 + t85;
t29 = -qJD(4) * t74 + (-rSges(5,2) + t101) * V_base(5) + (-t99 - t102) * t75 + t84;
t28 = t99 * V_base(4) + (t88 - t100) * V_base(5) + t86;
t1 = m(1) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(2) * (t37 ^ 2 + t38 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t34 ^ 2 + t35 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t31 ^ 2 + t32 ^ 2 + t33 ^ 2) / 0.2e1 + m(5) * (t28 ^ 2 + t29 ^ 2 + t30 ^ 2) / 0.2e1 + (Icges(1,2) / 0.2e1 + Icges(2,1) / 0.2e1 - t82 * t64 / 0.2e1 + t66 * t104 + t89 * t74 + t91 * t73) * V_base(5) ^ 2 + ((Icges(1,1) / 0.2e1 + Icges(2,3) / 0.2e1 + t63 * t104 + t82 * t105 + t92 * t74 + t90 * t73) * V_base(4) + (Icges(1,4) + Icges(2,5) + (t64 / 0.2e1 + t105) * t83 + (t66 / 0.2e1 - t63 / 0.2e1) * t82 + (t90 - t91) * t74 + (t89 - t92) * t73) * V_base(5)) * V_base(4) + ((Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(6) + (Icges(2,4) + Icges(1,6)) * V_base(5) + (Icges(1,5) + Icges(2,6)) * V_base(4)) * V_base(6) + ((Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t75 + (-Icges(3,5) * t83 + Icges(3,6) * t82) * V_base(5) + (-Icges(3,5) * t82 - Icges(3,6) * t83) * V_base(4) + (-t108 * V_base(4) + t109 * V_base(5)) * t74 + (t108 * V_base(5) + t109 * V_base(4)) * t73) * t75;
T  = t1;
