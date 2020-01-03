% Calculate kinetic energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:16
% EndTime: 2019-12-31 18:02:16
% DurationCPUTime: 0.18s
% Computational Cost: add. (171->59), mult. (309->93), div. (0->0), fcn. (130->6), ass. (0->27)
t88 = m(3) / 0.2e1;
t71 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t85 = qJ(2) * qJD(1);
t67 = t77 * t71 + t78 * t85;
t65 = -qJD(1) * pkin(6) + t67;
t81 = sin(qJ(4));
t83 = cos(qJ(4));
t62 = t81 * qJD(3) + t83 * t65;
t87 = qJD(1) * t81;
t86 = t83 * qJD(1);
t66 = t78 * t71 - t77 * t85;
t64 = qJD(1) * pkin(3) - t66;
t61 = t83 * qJD(3) - t81 * t65;
t82 = cos(qJ(5));
t80 = sin(qJ(5));
t75 = -qJD(1) * pkin(1) + qJD(2);
t72 = qJD(5) + t86;
t69 = t80 * qJD(4) - t82 * t87;
t68 = t82 * qJD(4) + t80 * t87;
t60 = (pkin(4) * t83 + pkin(7) * t81) * qJD(1) + t64;
t59 = qJD(4) * pkin(7) + t62;
t58 = -qJD(4) * pkin(4) - t61;
t57 = t82 * t59 + t80 * t60;
t56 = -t80 * t59 + t82 * t60;
t1 = m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t64 ^ 2) / 0.2e1 + t75 ^ 2 * t88 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t72 / 0.2e1) * t72 + (t61 * mrSges(5,1) - t62 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t58 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t72 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t58 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t72 + Ifges(6,2) * t68 / 0.2e1) * t68 + (-t75 * mrSges(3,1) - t66 * mrSges(4,1) + t67 * mrSges(4,2) + (t64 * mrSges(5,1) - t62 * mrSges(5,3) - Ifges(5,6) * qJD(4) + Ifges(5,2) * t86 / 0.2e1) * t83 + (-t64 * mrSges(5,2) + t61 * mrSges(5,3) - Ifges(5,5) * qJD(4)) * t81 + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (qJ(2) * t88 + mrSges(3,3)) * qJ(2) + (Ifges(5,4) * t83 + Ifges(5,1) * t81 / 0.2e1) * t81) * qJD(1)) * qJD(1);
T = t1;
