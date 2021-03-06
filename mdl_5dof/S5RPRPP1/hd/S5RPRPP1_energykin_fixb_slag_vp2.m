% Calculate kinetic energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:35
% EndTime: 2019-12-31 18:08:36
% DurationCPUTime: 0.30s
% Computational Cost: add. (184->70), mult. (428->103), div. (0->0), fcn. (230->6), ass. (0->25)
t84 = m(3) / 0.2e1;
t75 = sin(pkin(7));
t69 = (pkin(1) * t75 + pkin(6)) * qJD(1);
t78 = cos(qJ(3));
t73 = t78 * qJD(2);
t77 = sin(qJ(3));
t82 = qJ(4) * qJD(1);
t61 = qJD(3) * pkin(3) + t73 + (-t69 - t82) * t77;
t64 = t77 * qJD(2) + t78 * t69;
t62 = t78 * t82 + t64;
t74 = sin(pkin(8));
t83 = cos(pkin(8));
t57 = t74 * t61 + t83 * t62;
t76 = cos(pkin(7));
t81 = -pkin(1) * t76 - pkin(2);
t56 = t83 * t61 - t74 * t62;
t65 = qJD(4) + (-pkin(3) * t78 + t81) * qJD(1);
t70 = t81 * qJD(1);
t67 = (t74 * t78 + t83 * t77) * qJD(1);
t66 = (t74 * t77 - t78 * t83) * qJD(1);
t63 = -t77 * t69 + t73;
t58 = t66 * pkin(4) - t67 * qJ(5) + t65;
t55 = qJD(3) * qJ(5) + t57;
t54 = -qJD(3) * pkin(4) + qJD(5) - t56;
t1 = qJD(2) ^ 2 * t84 + m(4) * (t63 ^ 2 + t64 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t65 ^ 2) / 0.2e1 + (t65 * mrSges(5,2) + t54 * mrSges(6,2) - t56 * mrSges(5,3) - t58 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t67) * t67 + (t65 * mrSges(5,1) + t58 * mrSges(6,1) - t55 * mrSges(6,2) - t57 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t66 + (-Ifges(5,4) + Ifges(6,5)) * t67) * t66 + (t63 * mrSges(4,1) + t56 * mrSges(5,1) - t54 * mrSges(6,1) - t64 * mrSges(4,2) - t57 * mrSges(5,2) + t55 * mrSges(6,3) + (Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (Ifges(6,4) + Ifges(5,5)) * t67 + (-Ifges(5,6) + Ifges(6,6)) * t66) * qJD(3) + (t70 * (-mrSges(4,1) * t78 + mrSges(4,2) * t77) + (-t63 * t77 + t64 * t78) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t77 + Ifges(4,6) * t78) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t76 * mrSges(3,1) - t75 * mrSges(3,2) + (t75 ^ 2 + t76 ^ 2) * t84 * pkin(1)) * pkin(1) + Ifges(4,2) * t78 ^ 2 / 0.2e1 + (Ifges(4,4) * t78 + Ifges(4,1) * t77 / 0.2e1) * t77) * qJD(1)) * qJD(1);
T = t1;
