% Calculate kinetic energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:46
% EndTime: 2020-01-03 12:01:47
% DurationCPUTime: 0.18s
% Computational Cost: add. (239->60), mult. (368->101), div. (0->0), fcn. (186->8), ass. (0->29)
t85 = qJD(1) + qJD(2);
t92 = cos(qJ(4));
t98 = t85 * t92;
t93 = cos(qJ(2));
t97 = qJD(1) * pkin(1);
t78 = t85 * pkin(2) + t93 * t97;
t86 = sin(pkin(9));
t87 = cos(pkin(9));
t90 = sin(qJ(2));
t96 = t90 * t97;
t74 = t86 * t78 + t87 * t96;
t72 = t85 * pkin(7) + t74;
t89 = sin(qJ(4));
t68 = t89 * qJD(3) + t92 * t72;
t73 = t87 * t78 - t86 * t96;
t91 = cos(qJ(5));
t88 = sin(qJ(5));
t84 = qJD(4) + qJD(5);
t82 = t92 * qJD(3);
t76 = (t88 * t92 + t89 * t91) * t85;
t75 = (-t88 * t89 + t91 * t92) * t85;
t71 = -t85 * pkin(3) - t73;
t69 = (-pkin(4) * t92 - pkin(3)) * t85 - t73;
t67 = -t89 * t72 + t82;
t66 = pkin(8) * t98 + t68;
t65 = qJD(4) * pkin(4) + t82 + (-pkin(8) * t85 - t72) * t89;
t64 = t88 * t65 + t91 * t66;
t63 = t91 * t65 - t88 * t66;
t1 = m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t90 ^ 2 + t93 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t76 / 0.2e1) * t76 + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t76 + Ifges(6,6) * t84 + Ifges(6,2) * t75 / 0.2e1) * t75 + (t73 * mrSges(4,1) - t74 * mrSges(4,2) + (mrSges(3,1) * t93 - mrSges(3,2) * t90) * t97 + (-t71 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t98 / 0.2e1) * t92 + (t71 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t89 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t92 + Ifges(5,1) * t89 / 0.2e1) * t89) * t85) * t85;
T = t1;
