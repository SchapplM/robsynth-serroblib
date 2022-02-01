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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:24:48
% EndTime: 2022-01-23 09:24:49
% DurationCPUTime: 0.46s
% Computational Cost: add. (351->81), mult. (900->136), div. (0->0), fcn. (604->8), ass. (0->40)
t101 = cos(qJ(3));
t95 = sin(pkin(8));
t97 = cos(pkin(8));
t83 = qJD(2) + (-pkin(2) * t97 - pkin(6) * t95 - pkin(1)) * qJD(1);
t82 = t101 * t83;
t107 = t97 * qJD(1);
t88 = qJD(3) - t107;
t99 = sin(qJ(3));
t73 = t88 * pkin(3) + t82 + (-qJ(2) * t97 * t99 - qJ(4) * t101 * t95) * qJD(1);
t108 = t95 * qJD(1);
t105 = t99 * t108;
t106 = qJD(1) * qJ(2);
t104 = t97 * t106;
t78 = t101 * t104 + t99 * t83;
t76 = -qJ(4) * t105 + t78;
t94 = sin(pkin(9));
t96 = cos(pkin(9));
t68 = t94 * t73 + t96 * t76;
t102 = qJD(1) ^ 2;
t109 = t102 * qJ(2) ^ 2;
t84 = pkin(3) * t105 + t95 * t106 + qJD(4);
t67 = t96 * t73 - t94 * t76;
t100 = cos(qJ(5));
t98 = sin(qJ(5));
t93 = t97 ^ 2;
t92 = t95 ^ 2;
t91 = -qJD(1) * pkin(1) + qJD(2);
t89 = t92 * t109;
t86 = qJD(5) + t88;
t80 = (t101 * t96 - t94 * t99) * t108;
t79 = (-t101 * t94 - t96 * t99) * t108;
t77 = -t104 * t99 + t82;
t75 = -t79 * pkin(4) + t84;
t70 = t100 * t80 + t98 * t79;
t69 = t100 * t79 - t98 * t80;
t66 = t79 * pkin(7) + t68;
t65 = t88 * pkin(4) - t80 * pkin(7) + t67;
t64 = t100 * t66 + t98 * t65;
t63 = t100 * t65 - t98 * t66;
t1 = m(4) * (t77 ^ 2 + t78 ^ 2 + t89) / 0.2e1 + m(3) * (t109 * t93 + t91 ^ 2 + t89) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t75 ^ 2) / 0.2e1 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t84 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,1) * t80 / 0.2e1) * t80 + (-t84 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t80 + Ifges(5,2) * t79 / 0.2e1) * t79 + (t75 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t75 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t86 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t77 * mrSges(4,1) + t67 * mrSges(5,1) - t78 * mrSges(4,2) - t68 * mrSges(5,2) + Ifges(5,5) * t80 + Ifges(5,6) * t79 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t88) * t88 + ((-t91 * mrSges(3,1) + Ifges(3,2) * t107 / 0.2e1) * t97 + (t91 * mrSges(3,2) + (Ifges(3,4) * t97 + (Ifges(3,1) / 0.2e1 + (qJ(2) * mrSges(4,1) + Ifges(4,2) * t99 / 0.2e1) * t99 + (qJ(2) * mrSges(4,2) - Ifges(4,4) * t99 + Ifges(4,1) * t101 / 0.2e1) * t101) * t95) * qJD(1) + (-t77 * t101 - t78 * t99) * mrSges(4,3) + t88 * (Ifges(4,5) * t101 - Ifges(4,6) * t99)) * t95) * qJD(1) + (Ifges(2,3) / 0.2e1 + (t92 + t93) * qJ(2) * mrSges(3,3)) * t102;
T = t1;
