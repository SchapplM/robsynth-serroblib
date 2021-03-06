% Calculate kinetic energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:16
% EndTime: 2019-12-31 20:03:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (204->80), mult. (440->107), div. (0->0), fcn. (216->4), ass. (0->23)
t87 = pkin(6) * mrSges(3,3);
t79 = sin(qJ(2));
t85 = t79 * qJD(1);
t84 = pkin(6) * t85 + qJD(3);
t62 = -pkin(7) * t85 + (-pkin(2) - pkin(3)) * qJD(2) + t84;
t81 = cos(qJ(2));
t86 = qJD(1) * t81;
t69 = pkin(6) * t86 + qJD(2) * qJ(3);
t66 = -pkin(7) * t86 + t69;
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t57 = t78 * t62 + t80 * t66;
t67 = -qJD(1) * pkin(1) - pkin(2) * t86 - qJ(3) * t85;
t56 = t80 * t62 - t78 * t66;
t61 = pkin(3) * t86 - t67;
t75 = -qJD(2) + qJD(4);
t68 = -qJD(2) * pkin(2) + t84;
t65 = (-t78 * t81 + t79 * t80) * qJD(1);
t64 = (-t78 * t79 - t80 * t81) * qJD(1);
t58 = -t64 * pkin(4) + qJD(5) + t61;
t55 = t64 * qJ(5) + t57;
t54 = t75 * pkin(4) - t65 * qJ(5) + t56;
t1 = m(4) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t61 ^ 2) / 0.2e1 + (-t68 * mrSges(4,1) + t69 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + (t56 * mrSges(5,1) + t54 * mrSges(6,1) - t57 * mrSges(5,2) - t55 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t75) * t75 + (t61 * mrSges(5,2) + t58 * mrSges(6,2) - t56 * mrSges(5,3) - t54 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t65 + (Ifges(5,5) + Ifges(6,5)) * t75) * t65 + (-t61 * mrSges(5,1) - t58 * mrSges(6,1) + t57 * mrSges(5,3) + t55 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t64 + (Ifges(5,6) + Ifges(6,6)) * t75 + (Ifges(5,4) + Ifges(6,4)) * t65) * t64 + ((-t67 * mrSges(4,1) + t69 * mrSges(4,2) + (-pkin(6) * mrSges(3,2) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t81 + (t68 * mrSges(4,2) - t67 * mrSges(4,3) + (-pkin(6) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t79 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t79 ^ 2 + t81 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t87 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t81) * t81 + (-pkin(1) * mrSges(3,2) + (t87 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t79 + (Ifges(3,4) - Ifges(4,5)) * t81) * t79) * qJD(1)) * qJD(1);
T = t1;
