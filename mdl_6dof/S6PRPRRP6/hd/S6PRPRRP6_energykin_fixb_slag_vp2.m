% Calculate kinetic energy for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:03
% EndTime: 2019-03-08 20:18:03
% DurationCPUTime: 0.28s
% Computational Cost: add. (236->79), mult. (464->114), div. (0->0), fcn. (262->8), ass. (0->33)
t102 = sin(qJ(2));
t98 = sin(pkin(6));
t111 = qJD(1) * t98;
t108 = t102 * t111;
t90 = qJD(2) * qJ(3) + t108;
t113 = t90 ^ 2;
t112 = cos(qJ(5));
t100 = sin(qJ(5));
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t99 = cos(pkin(6));
t110 = qJD(1) * t99;
t104 = cos(qJ(2));
t107 = -t104 * t111 + qJD(3);
t86 = (-pkin(2) - pkin(8)) * qJD(2) + t107;
t82 = t101 * t86 + t103 * t110;
t80 = qJD(4) * pkin(9) + t82;
t84 = t108 + (pkin(4) * t101 - pkin(9) * t103 + qJ(3)) * qJD(2);
t76 = t100 * t84 + t112 * t80;
t109 = qJD(2) * t103;
t81 = -t101 * t110 + t103 * t86;
t79 = -qJD(4) * pkin(4) - t81;
t75 = -t100 * t80 + t112 * t84;
t106 = qJD(1) ^ 2;
t95 = t99 ^ 2 * t106;
t94 = t101 * qJD(2) + qJD(5);
t89 = t100 * qJD(4) + t109 * t112;
t88 = -qJD(4) * t112 + t100 * t109;
t87 = -qJD(2) * pkin(2) + t107;
t77 = t88 * pkin(5) - t89 * qJ(6) + t79;
t74 = t94 * qJ(6) + t76;
t73 = -t94 * pkin(5) + qJD(6) - t75;
t1 = m(3) * (t95 + (t102 ^ 2 + t104 ^ 2) * t98 ^ 2 * t106) / 0.2e1 + m(2) * t106 / 0.2e1 + m(4) * (t87 ^ 2 + t113 + t95) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t113) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + (t81 * mrSges(5,1) - t82 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t75 * mrSges(6,1) - t73 * mrSges(7,1) - t76 * mrSges(6,2) + t74 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t94) * t94 + (t79 * mrSges(6,2) + t73 * mrSges(7,2) - t75 * mrSges(6,3) - t77 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t89 + (Ifges(7,4) + Ifges(6,5)) * t94) * t89 + (t79 * mrSges(6,1) + t77 * mrSges(7,1) - t74 * mrSges(7,2) - t76 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t88 + (-Ifges(6,6) + Ifges(7,6)) * t94 + (-Ifges(6,4) + Ifges(7,5)) * t89) * t88 + (t87 * mrSges(4,2) + t90 * mrSges(4,3) + (mrSges(3,1) * t104 - mrSges(3,2) * t102) * t111 + (t90 * mrSges(5,2) - t81 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t109 / 0.2e1) * t103 + (t90 * mrSges(5,1) - t82 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t101 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t103 + Ifges(5,2) * t101 / 0.2e1) * t101) * qJD(2)) * qJD(2);
T  = t1;
