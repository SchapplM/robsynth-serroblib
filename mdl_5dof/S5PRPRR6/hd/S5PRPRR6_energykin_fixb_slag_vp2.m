% Calculate kinetic energy for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:22
% EndTime: 2019-12-05 15:56:22
% DurationCPUTime: 0.28s
% Computational Cost: add. (227->69), mult. (555->117), div. (0->0), fcn. (385->10), ass. (0->35)
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t92 = (t100 * t105 - t102 * t108) * qJD(2);
t115 = pkin(7) * qJD(2);
t106 = sin(qJ(2));
t101 = sin(pkin(5));
t114 = qJD(1) * t101;
t96 = qJD(2) * qJ(3) + t106 * t114;
t103 = cos(pkin(5));
t113 = qJD(1) * t103;
t98 = t102 * t113;
t84 = t98 + (-t96 - t115) * t100;
t89 = t100 * t113 + t102 * t96;
t85 = t102 * t115 + t89;
t80 = t105 * t84 + t108 * t85;
t109 = cos(qJ(2));
t112 = -t109 * t114 + qJD(3);
t79 = -t105 * t85 + t108 * t84;
t90 = (-pkin(3) * t102 - pkin(2)) * qJD(2) + t112;
t107 = cos(qJ(5));
t104 = sin(qJ(5));
t95 = -qJD(2) * pkin(2) + t112;
t93 = (t100 * t108 + t102 * t105) * qJD(2);
t91 = qJD(5) + t92;
t88 = -t100 * t96 + t98;
t87 = t104 * qJD(4) + t107 * t93;
t86 = t107 * qJD(4) - t104 * t93;
t81 = t92 * pkin(4) - t93 * pkin(8) + t90;
t78 = qJD(4) * pkin(8) + t80;
t77 = -qJD(4) * pkin(4) - t79;
t76 = t104 * t81 + t107 * t78;
t75 = -t104 * t78 + t107 * t81;
t1 = m(4) * (t88 ^ 2 + t89 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + (t90 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t91 / 0.2e1) * t91 - (-t90 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t93 - Ifges(5,2) * t92 / 0.2e1) * t92 + (t77 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t91 + Ifges(6,1) * t87 / 0.2e1) * t87 + (-t77 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t87 + Ifges(6,6) * t91 + Ifges(6,2) * t86 / 0.2e1) * t86 + (m(3) * (t103 ^ 2 + (t106 ^ 2 + t109 ^ 2) * t101 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,5) * t93 - Ifges(5,6) * t92 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t95 * (-mrSges(4,1) * t102 + mrSges(4,2) * t100) + (Ifges(4,2) * t102 ^ 2 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t102 + Ifges(4,1) * t100 / 0.2e1) * t100) * qJD(2) + (mrSges(3,1) * t109 - mrSges(3,2) * t106) * t114 + (-t88 * t100 + t89 * t102) * mrSges(4,3)) * qJD(2);
T = t1;
