% Calculate kinetic energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:52
% EndTime: 2019-03-09 03:30:53
% DurationCPUTime: 0.32s
% Computational Cost: add. (255->90), mult. (462->119), div. (0->0), fcn. (200->4), ass. (0->28)
t92 = cos(qJ(5));
t76 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t91 = t76 * mrSges(4,3);
t86 = cos(qJ(3));
t65 = qJD(4) + (pkin(4) * qJD(1) - t76) * t86 + (-pkin(3) - pkin(8)) * qJD(3);
t85 = sin(qJ(3));
t83 = qJD(1) * qJ(2);
t89 = qJD(1) * t85;
t90 = pkin(3) * t89 + t83;
t67 = (pkin(8) * t85 - qJ(4) * t86) * qJD(1) + t90;
t84 = sin(qJ(5));
t62 = t84 * t65 + t92 * t67;
t71 = -qJD(3) * qJ(4) - t85 * t76;
t88 = t86 * qJD(1);
t68 = -pkin(4) * t89 - t71;
t61 = t65 * t92 - t84 * t67;
t87 = qJD(1) ^ 2;
t82 = t87 * qJ(2) ^ 2;
t80 = -qJD(1) * pkin(1) + qJD(2);
t78 = qJD(5) + t88;
t73 = qJD(3) * t92 + t84 * t89;
t72 = t84 * qJD(3) - t89 * t92;
t70 = -qJ(4) * t88 + t90;
t69 = -qJD(3) * pkin(3) - t86 * t76 + qJD(4);
t63 = t72 * pkin(5) - t73 * qJ(6) + t68;
t60 = t78 * qJ(6) + t62;
t59 = -t78 * pkin(5) + qJD(6) - t61;
t1 = m(7) * (t59 ^ 2 + t60 ^ 2 + t63 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(3) * (t80 ^ 2 + t82) / 0.2e1 + m(4) * (t82 + (t85 ^ 2 + t86 ^ 2) * t76 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t87 + (t61 * mrSges(6,1) - t59 * mrSges(7,1) - t62 * mrSges(6,2) + t60 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t78) * t78 + (t69 * mrSges(5,2) - t71 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3) + (t86 * mrSges(4,1) - t85 * mrSges(4,2)) * t76) * qJD(3) + (t80 * mrSges(3,2) + (mrSges(4,2) * t83 + t69 * mrSges(5,1) - t70 * mrSges(5,3) + (-t91 + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(1)) * t86 + (-Ifges(5,4) + Ifges(4,5)) * qJD(3)) * t86 + (t71 * mrSges(5,1) - t70 * mrSges(5,2) + (qJ(2) * mrSges(4,1) + (-Ifges(4,4) - Ifges(5,6)) * t86) * qJD(1) + (-t91 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t85 + (Ifges(5,5) - Ifges(4,6)) * qJD(3)) * t85) * qJD(1) + (t68 * mrSges(6,2) + t59 * mrSges(7,2) - t61 * mrSges(6,3) - t63 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t73 + (Ifges(7,4) + Ifges(6,5)) * t78) * t73 + (t68 * mrSges(6,1) + t63 * mrSges(7,1) - t60 * mrSges(7,2) - t62 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t72 + (-Ifges(6,6) + Ifges(7,6)) * t78 + (-Ifges(6,4) + Ifges(7,5)) * t73) * t72;
T  = t1;
