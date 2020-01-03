% Calculate kinetic energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:39
% EndTime: 2019-12-31 19:05:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (225->57), mult. (310->93), div. (0->0), fcn. (130->6), ass. (0->27)
t89 = m(3) / 0.2e1;
t77 = -qJD(1) + qJD(3);
t88 = t77 / 0.2e1;
t72 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t86 = qJ(2) * qJD(1);
t70 = t80 * t72 + t83 * t86;
t66 = pkin(7) * t77 + t70;
t87 = t66 * mrSges(5,3);
t85 = pkin(8) * t77 + t66;
t69 = t72 * t83 - t80 * t86;
t82 = cos(qJ(4));
t81 = cos(qJ(5));
t79 = sin(qJ(4));
t78 = sin(qJ(5));
t76 = qJD(4) + qJD(5);
t75 = -qJD(1) * pkin(1) + qJD(2);
t68 = (t78 * t82 + t79 * t81) * t77;
t67 = (-t78 * t79 + t81 * t82) * t77;
t65 = -pkin(3) * t77 - t69;
t63 = (-pkin(4) * t82 - pkin(3)) * t77 - t69;
t62 = t85 * t82;
t61 = qJD(4) * pkin(4) - t85 * t79;
t60 = t61 * t78 + t62 * t81;
t59 = t61 * t81 - t62 * t78;
t1 = m(6) * (t59 ^ 2 + t60 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + (t79 ^ 2 + t82 ^ 2) * t66 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2) / 0.2e1 + (-qJD(1) * mrSges(3,1) + t75 * t89) * t75 + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (qJ(2) * t89 + mrSges(3,3)) * qJ(2)) * qJD(1) ^ 2 + (t59 * mrSges(6,1) - t60 * mrSges(6,2) + Ifges(6,3) * t76 / 0.2e1) * t76 + (Ifges(5,3) * qJD(4) / 0.2e1 + (-t79 * mrSges(5,1) - t82 * mrSges(5,2)) * t66) * qJD(4) + (t63 * mrSges(6,2) - t59 * mrSges(6,3) + Ifges(6,5) * t76 + Ifges(6,1) * t68 / 0.2e1) * t68 + (-t63 * mrSges(6,1) + t60 * mrSges(6,3) + Ifges(6,4) * t68 + Ifges(6,6) * t76 + Ifges(6,2) * t67 / 0.2e1) * t67 + (-t70 * mrSges(4,2) + t69 * mrSges(4,1) + Ifges(4,3) * t88 + (-t65 * mrSges(5,1) + Ifges(5,6) * qJD(4) + (Ifges(5,2) * t88 + t87) * t82) * t82 + (Ifges(5,4) * t82 * t77 + t65 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t88 + t87) * t79) * t79) * t77;
T = t1;
