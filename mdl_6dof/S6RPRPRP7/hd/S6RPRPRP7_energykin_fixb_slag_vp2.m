% Calculate kinetic energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:07
% EndTime: 2019-03-09 03:21:08
% DurationCPUTime: 0.39s
% Computational Cost: add. (397->90), mult. (802->131), div. (0->0), fcn. (486->6), ass. (0->33)
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t96 = sin(qJ(3));
t98 = cos(qJ(3));
t84 = (t93 * t98 + t94 * t96) * qJD(1);
t88 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t103 = t88 * mrSges(4,3);
t102 = qJD(1) / 0.2e1;
t101 = -qJ(4) * qJD(1) + t88;
t81 = qJD(3) * pkin(3) + t101 * t98;
t82 = t101 * t96;
t74 = t93 * t81 + t94 * t82;
t70 = qJD(3) * pkin(8) + t74;
t85 = (-t93 * t96 + t94 * t98) * qJD(1);
t92 = qJD(1) * qJ(2);
t86 = t96 * qJD(1) * pkin(3) + qJD(4) + t92;
t75 = t84 * pkin(4) - t85 * pkin(8) + t86;
t95 = sin(qJ(5));
t97 = cos(qJ(5));
t66 = t97 * t70 + t95 * t75;
t65 = -t95 * t70 + t97 * t75;
t73 = t94 * t81 - t93 * t82;
t69 = -qJD(3) * pkin(4) - t73;
t99 = qJD(1) ^ 2;
t91 = t99 * qJ(2) ^ 2;
t90 = -qJD(1) * pkin(1) + qJD(2);
t83 = qJD(5) + t84;
t77 = t95 * qJD(3) + t97 * t85;
t76 = t97 * qJD(3) - t95 * t85;
t67 = -t76 * pkin(5) + qJD(6) + t69;
t64 = t76 * qJ(6) + t66;
t63 = t83 * pkin(5) - t77 * qJ(6) + t65;
t1 = m(4) * (t91 + (t96 ^ 2 + t98 ^ 2) * t88 ^ 2) / 0.2e1 + m(3) * (t90 ^ 2 + t91) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t99 + (t86 * mrSges(5,2) - t73 * mrSges(5,3) + Ifges(5,1) * t85 / 0.2e1) * t85 - (-t86 * mrSges(5,1) + t74 * mrSges(5,3) + Ifges(5,4) * t85 - Ifges(5,2) * t84 / 0.2e1) * t84 + (t65 * mrSges(6,1) + t63 * mrSges(7,1) - t66 * mrSges(6,2) - t64 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t83) * t83 + (t90 * mrSges(3,2) + (mrSges(4,2) * t92 + (Ifges(4,1) * t102 - t103) * t98) * t98 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t98) * qJD(1) + (Ifges(4,2) * t102 - t103) * t96) * t96) * qJD(1) + (t73 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,5) * t85 - Ifges(5,6) * t84 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (t98 * mrSges(4,1) - t96 * mrSges(4,2)) * t88 + (Ifges(4,5) * t98 - Ifges(4,6) * t96) * qJD(1)) * qJD(3) + (t69 * mrSges(6,2) + t67 * mrSges(7,2) - t65 * mrSges(6,3) - t63 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t77 + (Ifges(6,5) + Ifges(7,5)) * t83) * t77 + (-t69 * mrSges(6,1) - t67 * mrSges(7,1) + t66 * mrSges(6,3) + t64 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t76 + (Ifges(6,6) + Ifges(7,6)) * t83 + (Ifges(6,4) + Ifges(7,4)) * t77) * t76;
T  = t1;
