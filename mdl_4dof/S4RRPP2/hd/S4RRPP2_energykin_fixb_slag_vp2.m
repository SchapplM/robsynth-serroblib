% Calculate kinetic energy for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RRPP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:23
% EndTime: 2018-11-14 13:52:23
% DurationCPUTime: 0.05s
% Computational Cost: add. (49->22), mult. (75->34), div. (0->0), fcn. (12->2), ass. (0->10)
t50 = pkin(1) * qJD(1);
t46 = cos(qJ(2));
t49 = -t46 * t50 + qJD(3);
t45 = sin(qJ(2));
t43 = qJD(1) + qJD(2);
t41 = t43 * qJ(3) + t45 * t50;
t40 = t41 ^ 2;
t39 = -t43 * pkin(2) + t49;
t38 = (-pkin(2) - pkin(3)) * t43 + t49;
t1 = m(4) * (t39 ^ 2 + t40) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t38 ^ 2 + t40) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t45 ^ 2 + t46 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (-t39 * mrSges(4,1) - t38 * mrSges(5,1) + (mrSges(5,2) + mrSges(4,3)) * t41 + (mrSges(3,1) * t46 - mrSges(3,2) * t45) * t50 + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t43) * t43;
T  = t1;
