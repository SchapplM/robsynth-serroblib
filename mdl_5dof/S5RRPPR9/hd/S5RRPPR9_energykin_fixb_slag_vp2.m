% Calculate kinetic energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:36
% EndTime: 2019-12-31 19:40:37
% DurationCPUTime: 0.27s
% Computational Cost: add. (174->81), mult. (362->107), div. (0->0), fcn. (144->4), ass. (0->26)
t87 = -pkin(2) - pkin(3);
t86 = pkin(6) * mrSges(3,3);
t78 = cos(qJ(2));
t85 = qJD(1) * t78;
t66 = pkin(6) * t85 + qJD(2) * qJ(3);
t76 = sin(qJ(2));
t84 = t76 * qJD(1);
t83 = pkin(6) * t84 + qJD(3);
t82 = qJ(4) * qJD(1);
t62 = -qJD(1) * pkin(1) - pkin(2) * t85 - qJ(3) * t84;
t58 = pkin(3) * t85 + qJD(4) - t62;
t61 = t78 * t82 - t66;
t81 = -t76 * t82 + t83;
t77 = cos(qJ(5));
t75 = sin(qJ(5));
t67 = qJD(5) + t84;
t65 = -qJD(2) * pkin(2) + t83;
t64 = -t75 * qJD(2) - t77 * t85;
t63 = -t77 * qJD(2) + t75 * t85;
t60 = qJD(2) * pkin(4) - t61;
t59 = t87 * qJD(2) + t81;
t57 = (-pkin(7) + t87) * qJD(2) + t81;
t56 = (pkin(4) * t76 + pkin(7) * t78) * qJD(1) + t58;
t55 = t75 * t56 + t77 * t57;
t54 = t77 * t56 - t75 * t57;
t1 = m(6) * (t54 ^ 2 + t55 ^ 2 + t60 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t61 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t54 * mrSges(6,1) - t55 * mrSges(6,2) + Ifges(6,3) * t67 / 0.2e1) * t67 + (t60 * mrSges(6,2) - t54 * mrSges(6,3) + Ifges(6,5) * t67 + Ifges(6,1) * t64 / 0.2e1) * t64 + (-t60 * mrSges(6,1) + t55 * mrSges(6,3) + Ifges(6,4) * t64 + Ifges(6,6) * t67 + Ifges(6,2) * t63 / 0.2e1) * t63 + (-t65 * mrSges(4,1) - t61 * mrSges(5,1) + t59 * mrSges(5,2) + t66 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t62 * mrSges(4,1) + t66 * mrSges(4,2) - t58 * mrSges(5,2) + t61 * mrSges(5,3) + (-pkin(6) * mrSges(3,2) + Ifges(5,5) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t78 + (t58 * mrSges(5,1) + t65 * mrSges(4,2) - t62 * mrSges(4,3) - t59 * mrSges(5,3) + (-pkin(6) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5) + Ifges(5,6)) * qJD(2)) * t76 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t76 ^ 2 + t78 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + t86) * t78) * t78 + (-pkin(1) * mrSges(3,2) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + t86) * t76 + (Ifges(3,4) + Ifges(5,4) - Ifges(4,5)) * t78) * t76) * qJD(1)) * qJD(1);
T = t1;
