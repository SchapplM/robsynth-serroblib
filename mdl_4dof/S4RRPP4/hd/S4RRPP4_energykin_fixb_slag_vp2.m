% Calculate kinetic energy for
% S4RRPP4
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:50
% EndTime: 2019-12-31 16:58:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (84->56), mult. (199->71), div. (0->0), fcn. (62->2), ass. (0->15)
t59 = pkin(2) + pkin(3);
t58 = pkin(5) * mrSges(3,3);
t51 = cos(qJ(2));
t57 = pkin(5) * qJD(1);
t46 = qJD(2) * qJ(3) + t51 * t57;
t50 = sin(qJ(2));
t56 = t50 * t57 + qJD(3);
t55 = qJ(4) * qJD(1);
t54 = qJ(3) * t50 + pkin(1);
t45 = -qJD(2) * pkin(2) + t56;
t44 = (-pkin(2) * t51 - t54) * qJD(1);
t43 = -t51 * t55 + t46;
t42 = -t59 * qJD(2) - t50 * t55 + t56;
t41 = qJD(4) + (t59 * t51 + t54) * qJD(1);
t1 = m(4) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(5) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + (-t45 * mrSges(4,1) - t42 * mrSges(5,1) + t43 * mrSges(5,2) + t46 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t44 * mrSges(4,1) + t41 * mrSges(5,1) + t46 * mrSges(4,2) - t43 * mrSges(5,3) + (-pkin(5) * mrSges(3,2) + Ifges(3,6) - Ifges(4,6) + Ifges(5,6)) * qJD(2)) * t51 + (t45 * mrSges(4,2) + t41 * mrSges(5,2) - t44 * mrSges(4,3) - t42 * mrSges(5,3) + (-pkin(5) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5) - Ifges(5,5)) * qJD(2)) * t50 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t50 ^ 2 + t51 ^ 2) * pkin(5) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t58 + Ifges(5,2) / 0.2e1) * t51) * t51 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + t58 + Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t50 + (Ifges(3,4) - Ifges(5,4) - Ifges(4,5)) * t51) * t50) * qJD(1)) * qJD(1);
T = t1;
