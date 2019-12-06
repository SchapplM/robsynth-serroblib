% Calculate kinetic energy for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:43
% EndTime: 2019-12-05 17:14:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (260->75), mult. (576->125), div. (0->0), fcn. (387->10), ass. (0->36)
t102 = sin(qJ(4));
t103 = sin(qJ(3));
t106 = cos(qJ(4));
t107 = cos(qJ(3));
t90 = (t102 * t103 - t106 * t107) * qJD(2);
t104 = sin(qJ(2));
t99 = sin(pkin(5));
t114 = qJD(1) * t99;
t93 = qJD(2) * pkin(7) + t104 * t114;
t100 = cos(pkin(5));
t113 = qJD(1) * t100;
t96 = t107 * t113;
t82 = qJD(3) * pkin(3) + t96 + (-pkin(8) * qJD(2) - t93) * t103;
t112 = qJD(2) * t107;
t87 = t103 * t113 + t107 * t93;
t83 = pkin(8) * t112 + t87;
t78 = t102 * t82 + t106 * t83;
t108 = cos(qJ(2));
t111 = t108 * t114;
t77 = -t102 * t83 + t106 * t82;
t88 = -t111 + (-pkin(3) * t107 - pkin(2)) * qJD(2);
t105 = cos(qJ(5));
t101 = sin(qJ(5));
t98 = qJD(3) + qJD(4);
t94 = -qJD(2) * pkin(2) - t111;
t91 = (t102 * t107 + t103 * t106) * qJD(2);
t89 = qJD(5) + t90;
t86 = -t103 * t93 + t96;
t85 = t101 * t98 + t105 * t91;
t84 = -t101 * t91 + t105 * t98;
t79 = t90 * pkin(4) - t91 * pkin(9) + t88;
t76 = t98 * pkin(9) + t78;
t75 = -t98 * pkin(4) - t77;
t74 = t101 * t79 + t105 * t76;
t73 = -t101 * t76 + t105 * t79;
t1 = m(5) * (t77 ^ 2 + t78 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + (t77 * mrSges(5,1) - t78 * mrSges(5,2) + Ifges(5,3) * t98 / 0.2e1) * t98 + (t73 * mrSges(6,1) - t74 * mrSges(6,2) + Ifges(6,3) * t89 / 0.2e1) * t89 + (t86 * mrSges(4,1) - t87 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t88 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,5) * t98 + Ifges(5,1) * t91 / 0.2e1) * t91 + (t75 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,5) * t89 + Ifges(6,1) * t85 / 0.2e1) * t85 - (-t88 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,6) * t98 - Ifges(5,2) * t90 / 0.2e1) * t90 + (-t75 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t85 + Ifges(6,6) * t89 + Ifges(6,2) * t84 / 0.2e1) * t84 + (m(3) * (t100 ^ 2 + (t104 ^ 2 + t108 ^ 2) * t99 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t108 - mrSges(3,2) * t104) * t114 + (-t94 * mrSges(4,1) + t87 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t112 / 0.2e1) * t107 + (t94 * mrSges(4,2) - t86 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t107 + Ifges(4,1) * t103 / 0.2e1) * qJD(2)) * t103) * qJD(2);
T = t1;
