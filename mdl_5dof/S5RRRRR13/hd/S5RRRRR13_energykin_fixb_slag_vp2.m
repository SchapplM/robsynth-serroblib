% Calculate kinetic energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:18
% EndTime: 2024-09-27 17:30:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (423->60), mult. (554->112), div. (0->0), fcn. (324->10), ass. (0->35)
t103 = sin(qJ(4));
t107 = cos(qJ(4));
t101 = cos(pkin(5));
t104 = sin(qJ(3));
t108 = cos(qJ(3));
t105 = sin(qJ(2));
t116 = pkin(1) * qJD(1);
t112 = t105 * t116;
t109 = cos(qJ(2));
t99 = qJD(1) + qJD(2);
t94 = t99 * pkin(2) + t109 * t116;
t89 = -t104 * t112 + t108 * t94;
t98 = qJD(3) + t99;
t86 = t98 * pkin(3) + t89;
t114 = t101 * t86;
t100 = sin(pkin(5));
t115 = t100 * t98;
t90 = t104 * t94 + t108 * t112;
t85 = pkin(9) * t115 + t90;
t80 = t103 * t114 + t107 * t85;
t95 = t101 * t98 + qJD(4);
t113 = pkin(10) * t115;
t106 = cos(qJ(5));
t102 = sin(qJ(5));
t93 = qJD(5) + t95;
t88 = (t102 * t107 + t103 * t106) * t115;
t87 = (-t102 * t103 + t106 * t107) * t115;
t83 = t107 * t114;
t81 = (-pkin(4) * t107 * t98 - t86) * t100;
t79 = -t103 * t85 + t83;
t78 = t107 * t113 + t80;
t77 = t95 * pkin(4) + t83 + (-t85 - t113) * t103;
t76 = t102 * t77 + t106 * t78;
t75 = -t102 * t78 + t106 * t77;
t1 = m(5) * (t100 ^ 2 * t86 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2) / 0.2e1 + (Ifges(3,3) * t99 / 0.2e1 + (mrSges(3,1) * t109 - mrSges(3,2) * t105) * t116) * t99 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,3) * t95 / 0.2e1) * t95 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t93 / 0.2e1) * t93 + (Ifges(2,3) / 0.2e1 + m(3) * (t105 ^ 2 + t109 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t81 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t93 + Ifges(6,1) * t88 / 0.2e1) * t88 + (-t81 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t88 + Ifges(6,6) * t93 + Ifges(6,2) * t87 / 0.2e1) * t87 + (-t90 * mrSges(4,2) + t89 * mrSges(4,1) + Ifges(4,3) * t98 / 0.2e1 + ((-t86 * (-mrSges(5,1) * t107 + mrSges(5,2) * t103) + (Ifges(5,2) * t107 ^ 2 / 0.2e1 + (Ifges(5,4) * t107 + Ifges(5,1) * t103 / 0.2e1) * t103) * t98) * t100 + (-t79 * t103 + t80 * t107) * mrSges(5,3) + t95 * (Ifges(5,5) * t103 + Ifges(5,6) * t107)) * t100) * t98;
T = t1;
