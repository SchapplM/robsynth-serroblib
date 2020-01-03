% Calculate kinetic energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:51
% EndTime: 2019-12-31 19:34:52
% DurationCPUTime: 0.32s
% Computational Cost: add. (272->83), mult. (648->119), div. (0->0), fcn. (396->6), ass. (0->32)
t92 = pkin(3) + pkin(7);
t91 = pkin(6) * mrSges(3,3);
t90 = pkin(6) + qJ(3);
t80 = sin(qJ(2));
t88 = qJD(1) * t80;
t74 = qJD(2) * pkin(2) - t90 * t88;
t82 = cos(qJ(2));
t87 = qJD(1) * t82;
t75 = t90 * t87;
t78 = sin(pkin(8));
t89 = cos(pkin(8));
t65 = t78 * t74 + t89 * t75;
t63 = -qJD(2) * qJ(4) - t65;
t64 = t89 * t74 - t78 * t75;
t86 = qJD(4) - t64;
t76 = qJD(3) + (-pkin(2) * t82 - pkin(1)) * qJD(1);
t71 = (t78 * t82 + t89 * t80) * qJD(1);
t85 = -qJ(4) * t71 + t76;
t81 = cos(qJ(5));
t79 = sin(qJ(5));
t70 = t78 * t88 - t89 * t87;
t69 = qJD(5) + t71;
t67 = qJD(2) * t81 + t70 * t79;
t66 = -qJD(2) * t79 + t70 * t81;
t62 = -qJD(2) * pkin(3) + t86;
t61 = pkin(3) * t70 + t85;
t60 = -pkin(4) * t70 - t63;
t59 = t71 * pkin(4) - t92 * qJD(2) + t86;
t58 = t92 * t70 + t85;
t57 = t58 * t81 + t59 * t79;
t56 = -t58 * t79 + t59 * t81;
t1 = m(4) * (t64 ^ 2 + t65 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t60 ^ 2) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t69 / 0.2e1) * t69 + (t60 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t69 + Ifges(6,1) * t67 / 0.2e1) * t67 + (-t60 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t67 + Ifges(6,6) * t69 + Ifges(6,2) * t66 / 0.2e1) * t66 + (t62 * mrSges(5,1) + t76 * mrSges(4,2) - t64 * mrSges(4,3) - t61 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t71) * t71 + (t76 * mrSges(4,1) + t63 * mrSges(5,1) - t61 * mrSges(5,2) - t65 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t70 + (-Ifges(4,4) - Ifges(5,6)) * t71) * t70 + (t64 * mrSges(4,1) - t65 * mrSges(4,2) + t62 * mrSges(5,2) - t63 * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (-Ifges(5,4) + Ifges(4,5)) * t71 + (Ifges(5,5) - Ifges(4,6)) * t70 + (Ifges(3,5) * t80 + Ifges(3,6) * t82 + (-mrSges(3,1) * t80 - mrSges(3,2) * t82) * pkin(6)) * qJD(1)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t80 ^ 2 + t82 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t91) * t82) * t82 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t82 + (Ifges(3,1) / 0.2e1 + t91) * t80) * t80) * qJD(1) ^ 2;
T = t1;
