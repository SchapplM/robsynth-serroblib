% Calculate kinetic energy for
% S5RRPRR15
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR15_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:05
% EndTime: 2019-12-31 20:41:05
% DurationCPUTime: 0.37s
% Computational Cost: add. (284->84), mult. (580->123), div. (0->0), fcn. (314->6), ass. (0->32)
t97 = -pkin(2) - pkin(7);
t96 = pkin(6) * mrSges(3,3);
t90 = cos(qJ(2));
t87 = sin(qJ(2));
t93 = -qJ(3) * t87 - pkin(1);
t70 = (t97 * t90 + t93) * qJD(1);
t83 = t87 * qJD(1);
t94 = pkin(6) * t83 + qJD(3);
t71 = pkin(3) * t83 + t97 * qJD(2) + t94;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t63 = t89 * t70 + t86 * t71;
t95 = t90 * qJD(1);
t77 = -pkin(6) * t95 - qJD(2) * qJ(3);
t79 = t83 + qJD(4);
t72 = pkin(3) * t95 - t77;
t62 = -t70 * t86 + t89 * t71;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t78 = qJD(5) + t79;
t76 = -qJD(2) * pkin(2) + t94;
t75 = qJD(2) * t89 - t86 * t95;
t74 = -qJD(2) * t86 - t89 * t95;
t73 = (-pkin(2) * t90 + t93) * qJD(1);
t66 = -pkin(4) * t74 + t72;
t65 = t74 * t85 + t75 * t88;
t64 = t74 * t88 - t75 * t85;
t61 = pkin(8) * t74 + t63;
t60 = pkin(4) * t79 - pkin(8) * t75 + t62;
t59 = t60 * t85 + t61 * t88;
t58 = t60 * t88 - t61 * t85;
t1 = m(4) * (t73 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t58 ^ 2 + t59 ^ 2 + t66 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t72 ^ 2) / 0.2e1 + (t62 * mrSges(5,1) - t63 * mrSges(5,2) + Ifges(5,3) * t79 / 0.2e1) * t79 + (t58 * mrSges(6,1) - t59 * mrSges(6,2) + Ifges(6,3) * t78 / 0.2e1) * t78 + (t72 * mrSges(5,2) - t62 * mrSges(5,3) + Ifges(5,5) * t79 + Ifges(5,1) * t75 / 0.2e1) * t75 + (t66 * mrSges(6,2) - t58 * mrSges(6,3) + Ifges(6,5) * t78 + Ifges(6,1) * t65 / 0.2e1) * t65 + (-t72 * mrSges(5,1) + t63 * mrSges(5,3) + Ifges(5,4) * t75 + Ifges(5,6) * t79 + Ifges(5,2) * t74 / 0.2e1) * t74 + (-t66 * mrSges(6,1) + t59 * mrSges(6,3) + Ifges(6,4) * t65 + Ifges(6,6) * t78 + Ifges(6,2) * t64 / 0.2e1) * t64 + (t76 * mrSges(4,2) - t77 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + ((-t77 * mrSges(4,1) + t73 * mrSges(4,2) + (-pkin(6) * mrSges(3,2) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t90 + (t76 * mrSges(4,1) - t73 * mrSges(4,3) + (-pkin(6) * mrSges(3,1) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t87 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t87 ^ 2 + t90 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t96) * t90) * t90 + (-pkin(1) * mrSges(3,2) + (Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1 + t96) * t87 + (Ifges(3,4) + Ifges(4,6)) * t90) * t87) * qJD(1)) * qJD(1);
T = t1;
