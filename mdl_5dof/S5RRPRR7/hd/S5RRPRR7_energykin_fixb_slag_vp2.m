% Calculate kinetic energy for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:26
% EndTime: 2019-12-31 20:15:26
% DurationCPUTime: 0.19s
% Computational Cost: add. (212->55), mult. (273->93), div. (0->0), fcn. (114->6), ass. (0->26)
t72 = qJD(1) + qJD(2);
t75 = sin(qJ(2));
t84 = pkin(1) * qJD(1);
t83 = t75 * t84;
t68 = qJ(3) * t72 + t83;
t87 = t68 ^ 2;
t86 = t72 / 0.2e1;
t78 = cos(qJ(2));
t81 = -t78 * t84 + qJD(3);
t66 = (-pkin(2) - pkin(7)) * t72 + t81;
t85 = t66 * mrSges(5,3);
t82 = -pkin(8) * t72 + t66;
t77 = cos(qJ(4));
t76 = cos(qJ(5));
t74 = sin(qJ(4));
t73 = sin(qJ(5));
t71 = qJD(4) + qJD(5);
t67 = -pkin(2) * t72 + t81;
t64 = t83 + (pkin(4) * t74 + qJ(3)) * t72;
t63 = (-t73 * t74 + t76 * t77) * t72;
t62 = (-t73 * t77 - t74 * t76) * t72;
t61 = t82 * t74;
t60 = qJD(4) * pkin(4) + t82 * t77;
t59 = t60 * t73 + t61 * t76;
t58 = t60 * t76 - t61 * t73;
t1 = m(6) * (t58 ^ 2 + t59 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t67 ^ 2 + t87) / 0.2e1 + m(5) * (t87 + (t74 ^ 2 + t77 ^ 2) * t66 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t75 ^ 2 + t78 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t58 * mrSges(6,1) - t59 * mrSges(6,2) + Ifges(6,3) * t71 / 0.2e1) * t71 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t77 * mrSges(5,1) - t74 * mrSges(5,2)) * t66) * qJD(4) + (t64 * mrSges(6,2) - t58 * mrSges(6,3) + Ifges(6,5) * t71 + Ifges(6,1) * t63 / 0.2e1) * t63 + (-t64 * mrSges(6,1) + t59 * mrSges(6,3) + Ifges(6,4) * t63 + Ifges(6,6) * t71 + Ifges(6,2) * t62 / 0.2e1) * t62 + (t67 * mrSges(4,2) + t68 * mrSges(4,3) + (mrSges(3,1) * t78 - mrSges(3,2) * t75) * t84 + (t68 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t86 - t85) * t77) * t77 + (t68 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t86 - t85) * t74) * t74 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(5,4) * t77 * t74) * t72) * t72;
T = t1;
