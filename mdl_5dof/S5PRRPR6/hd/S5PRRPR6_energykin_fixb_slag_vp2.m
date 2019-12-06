% Calculate kinetic energy for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:00
% EndTime: 2019-12-05 16:30:00
% DurationCPUTime: 0.28s
% Computational Cost: add. (252->74), mult. (578->121), div. (0->0), fcn. (383->10), ass. (0->35)
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t103 = cos(pkin(5));
t114 = qJD(1) * t103;
t106 = sin(qJ(2));
t101 = sin(pkin(5));
t115 = qJD(1) * t101;
t95 = qJD(2) * pkin(7) + t106 * t115;
t89 = t105 * t114 + t108 * t95;
t87 = qJD(3) * qJ(4) + t89;
t109 = cos(qJ(2));
t111 = t109 * t115;
t90 = -t111 + (-pkin(3) * t108 - qJ(4) * t105 - pkin(2)) * qJD(2);
t79 = t100 * t90 + t102 * t87;
t113 = qJD(2) * t105;
t112 = t108 * qJD(2);
t78 = -t100 * t87 + t102 * t90;
t88 = -t105 * t95 + t108 * t114;
t84 = -qJD(3) * pkin(3) + qJD(4) - t88;
t107 = cos(qJ(5));
t104 = sin(qJ(5));
t98 = qJD(5) - t112;
t96 = -qJD(2) * pkin(2) - t111;
t94 = t100 * qJD(3) + t102 * t113;
t93 = t102 * qJD(3) - t100 * t113;
t82 = t104 * t93 + t107 * t94;
t81 = -t104 * t94 + t107 * t93;
t80 = -t93 * pkin(4) + t84;
t77 = t93 * pkin(8) + t79;
t76 = -pkin(4) * t112 - t94 * pkin(8) + t78;
t75 = t104 * t76 + t107 * t77;
t74 = -t104 * t77 + t107 * t76;
t1 = m(4) * (t88 ^ 2 + t89 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t80 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t84 ^ 2) / 0.2e1 + (t74 * mrSges(6,1) - t75 * mrSges(6,2) + Ifges(6,3) * t98 / 0.2e1) * t98 + (t84 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t94 / 0.2e1) * t94 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t84 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,2) * t93 / 0.2e1) * t93 + (t80 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,5) * t98 + Ifges(6,1) * t82 / 0.2e1) * t82 + (-t80 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t82 + Ifges(6,6) * t98 + Ifges(6,2) * t81 / 0.2e1) * t81 + (m(3) * (t103 ^ 2 + (t106 ^ 2 + t109 ^ 2) * t101 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t109 - mrSges(3,2) * t106) * t115 + (t96 * mrSges(4,2) - t88 * mrSges(4,3) + Ifges(4,5) * qJD(3) + Ifges(4,1) * t113 / 0.2e1) * t105 + (-t96 * mrSges(4,1) - t78 * mrSges(5,1) + t79 * mrSges(5,2) + t89 * mrSges(4,3) - Ifges(5,5) * t94 + Ifges(4,6) * qJD(3) - Ifges(5,6) * t93 + (Ifges(4,4) * t105 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t108) * qJD(2)) * t108) * qJD(2);
T = t1;
