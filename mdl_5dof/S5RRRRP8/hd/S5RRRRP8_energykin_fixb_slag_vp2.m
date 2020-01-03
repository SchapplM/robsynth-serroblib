% Calculate kinetic energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:35
% EndTime: 2019-12-31 21:59:35
% DurationCPUTime: 0.38s
% Computational Cost: add. (368->83), mult. (780->121), div. (0->0), fcn. (502->6), ass. (0->29)
t93 = pkin(6) * mrSges(3,3);
t85 = sin(qJ(2));
t88 = cos(qJ(2));
t74 = (-pkin(2) * t88 - pkin(7) * t85 - pkin(1)) * qJD(1);
t91 = t88 * qJD(1);
t79 = pkin(6) * t91 + qJD(2) * pkin(7);
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t69 = t87 * t74 - t79 * t84;
t92 = t85 * qJD(1);
t76 = qJD(2) * t84 + t87 * t92;
t81 = qJD(3) - t91;
t64 = pkin(3) * t81 - pkin(8) * t76 + t69;
t70 = t84 * t74 + t87 * t79;
t75 = qJD(2) * t87 - t84 * t92;
t66 = pkin(8) * t75 + t70;
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t60 = t83 * t64 + t86 * t66;
t59 = t86 * t64 - t66 * t83;
t78 = -qJD(2) * pkin(2) + pkin(6) * t92;
t71 = -pkin(3) * t75 + t78;
t80 = qJD(4) + t81;
t68 = t75 * t83 + t76 * t86;
t67 = t75 * t86 - t76 * t83;
t61 = -pkin(4) * t67 + qJD(5) + t71;
t58 = qJ(5) * t67 + t60;
t57 = pkin(4) * t80 - qJ(5) * t68 + t59;
t1 = m(6) * (t57 ^ 2 + t58 ^ 2 + t61 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t71 ^ 2) / 0.2e1 + (t69 * mrSges(4,1) - t70 * mrSges(4,2) + Ifges(4,3) * t81 / 0.2e1) * t81 + (t78 * mrSges(4,2) - t69 * mrSges(4,3) + Ifges(4,5) * t81 + Ifges(4,1) * t76 / 0.2e1) * t76 + (-t78 * mrSges(4,1) + t70 * mrSges(4,3) + Ifges(4,4) * t76 + Ifges(4,6) * t81 + Ifges(4,2) * t75 / 0.2e1) * t75 + (t59 * mrSges(5,1) + t57 * mrSges(6,1) - t60 * mrSges(5,2) - t58 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t80) * t80 + (t71 * mrSges(5,2) + t61 * mrSges(6,2) - t59 * mrSges(5,3) - t57 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t68 + (Ifges(5,5) + Ifges(6,5)) * t80) * t68 + (-t71 * mrSges(5,1) - t61 * mrSges(6,1) + t60 * mrSges(5,3) + t58 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t67 + (Ifges(5,6) + Ifges(6,6)) * t80 + (Ifges(5,4) + Ifges(6,4)) * t68) * t67 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t85 ^ 2 + t88 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t93 + Ifges(3,2) / 0.2e1) * t88) * t88 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t88 + (t93 + Ifges(3,1) / 0.2e1) * t85) * t85) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t88 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t85) * qJD(2)) * qJD(1);
T = t1;
