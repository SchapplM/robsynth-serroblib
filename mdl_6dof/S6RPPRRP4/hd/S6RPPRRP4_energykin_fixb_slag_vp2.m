% Calculate kinetic energy for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:05
% EndTime: 2019-03-09 02:05:06
% DurationCPUTime: 0.28s
% Computational Cost: add. (289->79), mult. (506->111), div. (0->0), fcn. (236->6), ass. (0->30)
t101 = m(3) / 0.2e1;
t84 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t90 = sin(pkin(9));
t91 = cos(pkin(9));
t98 = qJ(2) * qJD(1);
t80 = t90 * t84 + t91 * t98;
t78 = -qJD(1) * pkin(7) + t80;
t94 = sin(qJ(4));
t96 = cos(qJ(4));
t74 = t94 * qJD(3) + t96 * t78;
t71 = qJD(4) * pkin(8) + t74;
t79 = t84 * t91 - t90 * t98;
t77 = qJD(1) * pkin(3) - t79;
t72 = (pkin(4) * t96 + pkin(8) * t94) * qJD(1) + t77;
t93 = sin(qJ(5));
t95 = cos(qJ(5));
t66 = t95 * t71 + t93 * t72;
t100 = t94 * qJD(1);
t99 = t96 * qJD(1);
t73 = qJD(3) * t96 - t94 * t78;
t65 = -t71 * t93 + t72 * t95;
t70 = -qJD(4) * pkin(4) - t73;
t88 = -qJD(1) * pkin(1) + qJD(2);
t85 = qJD(5) + t99;
t82 = qJD(4) * t93 - t100 * t95;
t81 = qJD(4) * t95 + t100 * t93;
t67 = -pkin(5) * t81 - qJ(6) * t82 + t70;
t64 = qJ(6) * t85 + t66;
t63 = -pkin(5) * t85 + qJD(6) - t65;
t1 = t88 ^ 2 * t101 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + (t73 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t65 * mrSges(6,1) - t63 * mrSges(7,1) - t66 * mrSges(6,2) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t85) * t85 + (t70 * mrSges(6,2) + t63 * mrSges(7,2) - t65 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t82 + (Ifges(7,4) + Ifges(6,5)) * t85) * t82 + (-t70 * mrSges(6,1) - t67 * mrSges(7,1) + t64 * mrSges(7,2) + t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t81 + (Ifges(6,6) - Ifges(7,6)) * t85 + (Ifges(6,4) - Ifges(7,5)) * t82) * t81 + (-t88 * mrSges(3,1) - t79 * mrSges(4,1) + t80 * mrSges(4,2) + (t77 * mrSges(5,1) - t74 * mrSges(5,3) - Ifges(5,6) * qJD(4) + Ifges(5,2) * t99 / 0.2e1) * t96 + (-t77 * mrSges(5,2) + t73 * mrSges(5,3) - Ifges(5,5) * qJD(4)) * t94 + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (qJ(2) * t101 + mrSges(3,3)) * qJ(2) + (Ifges(5,4) * t96 + Ifges(5,1) * t94 / 0.2e1) * t94) * qJD(1)) * qJD(1);
T  = t1;
