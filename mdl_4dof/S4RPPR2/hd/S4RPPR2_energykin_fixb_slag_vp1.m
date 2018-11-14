% Calculate kinetic energy for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:25
% EndTime: 2018-11-14 13:47:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (70->49), mult. (116->72), div. (0->0), fcn. (72->6), ass. (0->30)
t173 = sin(pkin(6));
t175 = sin(qJ(1));
t181 = t175 * t173;
t176 = cos(qJ(1));
t180 = t176 * t173;
t162 = qJD(1) * (t176 * pkin(1) + t175 * qJ(2));
t179 = -qJD(2) * t176 + t162;
t177 = qJD(3) ^ 2;
t174 = cos(pkin(6));
t172 = qJD(1) - qJD(4);
t171 = pkin(6) + qJ(4);
t170 = qJD(2) * t175;
t169 = cos(t171);
t168 = sin(t171);
t167 = qJD(1) * t176 * pkin(2);
t166 = t174 * pkin(3) + pkin(2);
t165 = t176 * rSges(2,1) - t175 * rSges(2,2);
t164 = t175 * rSges(2,1) + t176 * rSges(2,2);
t163 = t175 * pkin(1) - t176 * qJ(2);
t161 = t175 * t174 - t180;
t160 = -t176 * t174 - t181;
t159 = -t176 * t168 + t175 * t169;
t158 = -t175 * t168 - t176 * t169;
t157 = qJD(1) * (t176 * rSges(3,1) + t175 * rSges(3,3)) + t179;
t156 = t170 + (-t175 * rSges(3,1) + t176 * rSges(3,3) - t163) * qJD(1);
t155 = t167 + qJD(1) * (-t160 * rSges(4,1) + t161 * rSges(4,2)) + t179;
t154 = t170 + (-t161 * rSges(4,1) - t160 * rSges(4,2) - t175 * pkin(2) - t163) * qJD(1);
t153 = t162 + t167 + qJD(1) * pkin(3) * t181 + t172 * (-t158 * rSges(5,1) + t159 * rSges(5,2)) + (-qJD(2) + qJD(1) * (-pkin(2) + t166)) * t176;
t152 = t170 - t172 * (t159 * rSges(5,1) + t158 * rSges(5,2)) + (pkin(3) * t180 - t175 * t166 - t163) * qJD(1);
t1 = m(3) * (t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(4) * (t154 ^ 2 + t155 ^ 2 + t177) / 0.2e1 + m(5) * (t152 ^ 2 + t153 ^ 2 + t177) / 0.2e1 + t172 ^ 2 * Icges(5,3) / 0.2e1 + (m(2) * (t164 ^ 2 + t165 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
