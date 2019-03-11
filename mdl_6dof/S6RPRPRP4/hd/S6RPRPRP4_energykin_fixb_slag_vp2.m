% Calculate kinetic energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:03
% EndTime: 2019-03-09 03:11:03
% DurationCPUTime: 0.40s
% Computational Cost: add. (266->90), mult. (531->122), div. (0->0), fcn. (256->6), ass. (0->32)
t102 = m(3) / 0.2e1;
t101 = -pkin(3) - pkin(8);
t88 = sin(pkin(9));
t82 = (pkin(1) * t88 + pkin(7)) * qJD(1);
t91 = sin(qJ(3));
t93 = cos(qJ(3));
t76 = t93 * qJD(2) - t91 * t82;
t97 = qJD(4) - t76;
t99 = t91 * qJD(1);
t69 = pkin(4) * t99 + t101 * qJD(3) + t97;
t89 = cos(pkin(9));
t98 = -pkin(1) * t89 - pkin(2);
t96 = -qJ(4) * t91 + t98;
t72 = (t101 * t93 + t96) * qJD(1);
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t66 = t90 * t69 + t92 * t72;
t77 = t91 * qJD(2) + t93 * t82;
t100 = qJD(1) * t93;
t74 = -qJD(3) * qJ(4) - t77;
t71 = pkin(4) * t100 - t74;
t65 = t92 * t69 - t90 * t72;
t84 = qJD(5) + t99;
t83 = t98 * qJD(1);
t81 = t92 * qJD(3) - t100 * t90;
t80 = t90 * qJD(3) + t100 * t92;
t75 = (-pkin(3) * t93 + t96) * qJD(1);
t73 = -qJD(3) * pkin(3) + t97;
t67 = t80 * pkin(5) - t81 * qJ(6) + t71;
t64 = t84 * qJ(6) + t66;
t63 = -t84 * pkin(5) + qJD(6) - t65;
t1 = qJD(2) ^ 2 * t102 + m(4) * (t76 ^ 2 + t77 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + (t65 * mrSges(6,1) - t63 * mrSges(7,1) - t66 * mrSges(6,2) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t84) * t84 + (t71 * mrSges(6,2) + t63 * mrSges(7,2) - t65 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81 + (Ifges(7,4) + Ifges(6,5)) * t84) * t81 + (t71 * mrSges(6,1) + t67 * mrSges(7,1) - t64 * mrSges(7,2) - t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t84 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + (t76 * mrSges(4,1) - t77 * mrSges(4,2) + t73 * mrSges(5,2) - t74 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3)) * qJD(3) + ((-t83 * mrSges(4,1) - t74 * mrSges(5,1) + t75 * mrSges(5,2) + t77 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t100) * t93 + (t73 * mrSges(5,1) + t83 * mrSges(4,2) - t76 * mrSges(4,3) - t75 * mrSges(5,3)) * t91 + ((-Ifges(5,5) + Ifges(4,6)) * t93 + (-Ifges(5,4) + Ifges(4,5)) * t91) * qJD(3) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t89 * mrSges(3,1) - t88 * mrSges(3,2) + (t88 ^ 2 + t89 ^ 2) * t102 * pkin(1)) * pkin(1) + ((Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t91 + (Ifges(4,4) + Ifges(5,6)) * t93) * t91) * qJD(1)) * qJD(1);
T  = t1;
