% Calculate kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:20
% EndTime: 2020-01-03 11:41:20
% DurationCPUTime: 0.51s
% Computational Cost: add. (351->81), mult. (900->136), div. (0->0), fcn. (604->8), ass. (0->40)
t92 = sin(pkin(8));
t94 = cos(pkin(8));
t80 = qJD(2) + (-pkin(2) * t94 - pkin(6) * t92 - pkin(1)) * qJD(1);
t98 = cos(qJ(3));
t79 = t98 * t80;
t104 = t94 * qJD(1);
t85 = qJD(3) - t104;
t96 = sin(qJ(3));
t70 = t85 * pkin(3) + t79 + (-qJ(2) * t94 * t96 - qJ(4) * t92 * t98) * qJD(1);
t105 = t92 * qJD(1);
t102 = t96 * t105;
t103 = qJD(1) * qJ(2);
t101 = t94 * t103;
t75 = t98 * t101 + t96 * t80;
t73 = -qJ(4) * t102 + t75;
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t65 = t91 * t70 + t93 * t73;
t99 = qJD(1) ^ 2;
t106 = qJ(2) ^ 2 * t99;
t81 = pkin(3) * t102 + t92 * t103 + qJD(4);
t64 = t93 * t70 - t91 * t73;
t97 = cos(qJ(5));
t95 = sin(qJ(5));
t90 = t94 ^ 2;
t89 = t92 ^ 2;
t88 = -qJD(1) * pkin(1) + qJD(2);
t86 = t89 * t106;
t83 = qJD(5) + t85;
t77 = (-t91 * t96 + t93 * t98) * t105;
t76 = (-t91 * t98 - t93 * t96) * t105;
t74 = -t101 * t96 + t79;
t72 = -t76 * pkin(4) + t81;
t67 = t95 * t76 + t97 * t77;
t66 = t97 * t76 - t95 * t77;
t63 = t76 * pkin(7) + t65;
t62 = t85 * pkin(4) - t77 * pkin(7) + t64;
t61 = t95 * t62 + t97 * t63;
t60 = t97 * t62 - t95 * t63;
t1 = m(6) * (t60 ^ 2 + t61 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t86) / 0.2e1 + m(3) * (t106 * t90 + t88 ^ 2 + t86) / 0.2e1 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t83 / 0.2e1) * t83 + (t81 * mrSges(5,2) - t64 * mrSges(5,3) + Ifges(5,1) * t77 / 0.2e1) * t77 + (-t81 * mrSges(5,1) + t65 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,2) * t76 / 0.2e1) * t76 + (t72 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t83 + Ifges(6,1) * t67 / 0.2e1) * t67 + (-t72 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t67 + Ifges(6,6) * t83 + Ifges(6,2) * t66 / 0.2e1) * t66 + (t74 * mrSges(4,1) + t64 * mrSges(5,1) - t75 * mrSges(4,2) - t65 * mrSges(5,2) + Ifges(5,5) * t77 + Ifges(5,6) * t76 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t85) * t85 + ((-t88 * mrSges(3,1) + Ifges(3,2) * t104 / 0.2e1) * t94 + (t88 * mrSges(3,2) + (Ifges(3,4) * t94 + (Ifges(3,1) / 0.2e1 + (qJ(2) * mrSges(4,2) + Ifges(4,1) * t98 / 0.2e1) * t98 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t98 + Ifges(4,2) * t96 / 0.2e1) * t96) * t92) * qJD(1) + (-t74 * t98 - t75 * t96) * mrSges(4,3) + t85 * (Ifges(4,5) * t98 - Ifges(4,6) * t96)) * t92) * qJD(1) + (Ifges(2,3) / 0.2e1 + (t89 + t90) * qJ(2) * mrSges(3,3)) * t99;
T = t1;
