% Calculate kinetic energy for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:54
% EndTime: 2019-12-31 16:21:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->28), mult. (86->44), div. (0->0), fcn. (20->2), ass. (0->10)
t51 = qJD(1) ^ 2;
t50 = qJD(2) ^ 2;
t49 = cos(qJ(4));
t48 = sin(qJ(4));
t47 = t50 * qJ(3) ^ 2;
t46 = -qJD(2) * pkin(2) + qJD(3);
t45 = qJD(3) + (-pkin(2) - pkin(5)) * qJD(2);
t44 = t49 * qJD(1) + t48 * t45;
t43 = -t48 * qJD(1) + t49 * t45;
t1 = m(5) * (t43 ^ 2 + t44 ^ 2 + t47) / 0.2e1 + m(4) * (t46 ^ 2 + t47 + t51) / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + qJ(3) * mrSges(4,3)) * t50 + (t43 * mrSges(5,1) - t44 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t46 * mrSges(4,2) + (-t43 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (qJ(3) * mrSges(5,2) + Ifges(5,1) * t49 / 0.2e1) * qJD(2)) * t49 + (-t44 * mrSges(5,3) - Ifges(5,6) * qJD(4) + (qJ(3) * mrSges(5,1) - Ifges(5,4) * t49 + Ifges(5,2) * t48 / 0.2e1) * qJD(2)) * t48) * qJD(2) + (m(2) + m(3)) * t51 / 0.2e1;
T = t1;
