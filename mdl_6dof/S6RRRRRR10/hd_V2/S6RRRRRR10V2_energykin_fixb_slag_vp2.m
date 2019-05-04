% Calculate kinetic energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10V2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:27
% EndTime: 2019-04-11 14:41:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (668->94), mult. (1311->159), div. (0->0), fcn. (1004->10), ass. (0->38)
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t106 = cos(qJ(3));
t107 = cos(qJ(2));
t90 = (t101 * t102 - t106 * t107) * qJD(1);
t104 = cos(qJ(5));
t100 = sin(qJ(4));
t105 = cos(qJ(4));
t91 = (t101 * t107 + t102 * t106) * qJD(1);
t96 = (-pkin(2) * t107 - pkin(1)) * qJD(1);
t85 = pkin(3) * t90 - pkin(5) * t91 + t96;
t112 = pkin(2) * qJD(2);
t97 = qJD(2) + qJD(3);
t94 = pkin(5) * t97 + t101 * t112;
t80 = t100 * t85 + t105 * t94;
t95 = -pkin(3) * t97 - t106 * t112;
t99 = sin(qJ(5));
t74 = -t104 * t95 + t80 * t99;
t115 = t74 ^ 2;
t78 = t100 * t94 - t105 * t85;
t114 = t78 ^ 2;
t76 = t104 * t80 + t99 * t95;
t88 = t100 * t97 + t105 * t91;
t89 = qJD(4) + t90;
t82 = t104 * t89 - t88 * t99;
t87 = -t100 * t91 + t105 * t97;
t103 = cos(qJ(6));
t98 = sin(qJ(6));
t86 = qJD(5) - t87;
t83 = t104 * t88 + t89 * t99;
t81 = qJD(6) - t82;
t73 = t103 * t83 + t86 * t98;
t72 = t103 * t86 - t83 * t98;
t71 = pkin(6) * t86 + t76;
t70 = -pkin(6) * t83 + t78;
t69 = t103 * t71 + t70 * t98;
t68 = t103 * t70 - t71 * t98;
t1 = m(7) * (t68 ^ 2 + t69 ^ 2 + t115) / 0.2e1 + m(6) * (t76 ^ 2 + t114 + t115) / 0.2e1 + Ifges(4,3) * t97 ^ 2 / 0.2e1 + m(4) * (t96 ^ 2 + (t101 ^ 2 + t106 ^ 2) * pkin(2) ^ 2 * qJD(2) ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t95 ^ 2 + t114) / 0.2e1 + (-t78 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,3) * t89 / 0.2e1) * t89 + (-t74 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t68 * mrSges(7,1) - t69 * mrSges(7,2) + Ifges(7,3) * t81 / 0.2e1) * t81 + (t96 * mrSges(4,2) + Ifges(4,5) * t97 + Ifges(4,1) * t91 / 0.2e1) * t91 + (t95 * mrSges(5,2) + t78 * mrSges(5,3) + Ifges(5,5) * t89 + Ifges(5,1) * t88 / 0.2e1) * t88 + (t78 * mrSges(6,2) + t74 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t83 / 0.2e1) * t83 + (t74 * mrSges(7,2) - t68 * mrSges(7,3) + Ifges(7,5) * t81 + Ifges(7,1) * t73 / 0.2e1) * t73 - (-t96 * mrSges(4,1) + Ifges(4,4) * t91 + Ifges(4,6) * t97 - Ifges(4,2) * t90 / 0.2e1) * t90 + (-t95 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t88 + Ifges(5,6) * t89 + Ifges(5,2) * t87 / 0.2e1) * t87 + (-t78 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t83 + Ifges(6,6) * t86 + Ifges(6,2) * t82 / 0.2e1) * t82 + (-t74 * mrSges(7,1) + t69 * mrSges(7,3) + Ifges(7,4) * t73 + Ifges(7,6) * t81 + Ifges(7,2) * t72 / 0.2e1) * t72 + (Ifges(3,3) * qJD(2) / 0.2e1 + (Ifges(3,5) * t102 + Ifges(3,6) * t107) * qJD(1) + (t101 * (-mrSges(4,2) * t97 - mrSges(4,3) * t90) + t106 * (mrSges(4,1) * t97 - mrSges(4,3) * t91)) * pkin(2)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * pkin(1) ^ 2 / 0.2e1 + (pkin(1) * mrSges(3,1) + Ifges(3,2) * t107 / 0.2e1) * t107 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t107 + Ifges(3,1) * t102 / 0.2e1) * t102) * qJD(1) ^ 2;
T  = t1;
