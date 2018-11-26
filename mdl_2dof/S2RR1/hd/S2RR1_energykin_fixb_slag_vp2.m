% Calculate kinetic energy for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S2RR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_energykin_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_energykin_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_energykin_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_energykin_fixb_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_energykin_fixb_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_energykin_fixb_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:19
% EndTime: 2018-11-16 16:44:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->11), mult. (50->24), div. (0->0), fcn. (14->2), ass. (0->4)
t21 = cos(qJ(2));
t26 = t21 ^ 2;
t20 = sin(qJ(2));
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + ((t26 * Ifges(3,2) / 0.2e1 + (Ifges(3,4) * t21 + Ifges(3,1) * t20 / 0.2e1) * t20 + Ifges(2,3) / 0.2e1 + (mrSges(3,3) + m(3) * pkin(1) / 0.2e1) * (t20 ^ 2 + t26) * pkin(1)) * qJD(1) + (-Ifges(3,5) * t20 - Ifges(3,6) * t21 + (t20 * mrSges(3,1) + t21 * mrSges(3,2)) * pkin(1)) * qJD(2)) * qJD(1);
T  = t1;
