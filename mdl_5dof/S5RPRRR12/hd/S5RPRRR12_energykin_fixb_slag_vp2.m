% Calculate kinetic energy for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:25
% EndTime: 2019-12-31 19:12:26
% DurationCPUTime: 0.32s
% Computational Cost: add. (253->70), mult. (481->114), div. (0->0), fcn. (266->6), ass. (0->31)
t82 = sin(qJ(4));
t83 = sin(qJ(3));
t85 = cos(qJ(4));
t86 = cos(qJ(3));
t71 = (t82 * t86 + t83 * t85) * qJD(1);
t75 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t91 = t75 * mrSges(4,3);
t90 = qJD(1) / 0.2e1;
t89 = -pkin(7) * qJD(1) + t75;
t68 = qJD(3) * pkin(3) + t89 * t86;
t69 = t89 * t83;
t62 = t82 * t68 + t85 * t69;
t80 = qJD(1) * qJ(2);
t73 = t83 * qJD(1) * pkin(3) + t80;
t61 = t68 * t85 - t69 * t82;
t87 = qJD(1) ^ 2;
t84 = cos(qJ(5));
t81 = sin(qJ(5));
t79 = t87 * qJ(2) ^ 2;
t78 = qJD(3) + qJD(4);
t77 = -qJD(1) * pkin(1) + qJD(2);
t72 = (-t82 * t83 + t85 * t86) * qJD(1);
t70 = qJD(5) + t71;
t65 = t72 * t84 + t78 * t81;
t64 = -t72 * t81 + t78 * t84;
t63 = pkin(4) * t71 - pkin(8) * t72 + t73;
t60 = pkin(8) * t78 + t62;
t59 = -pkin(4) * t78 - t61;
t58 = t60 * t84 + t63 * t81;
t57 = -t60 * t81 + t63 * t84;
t1 = m(4) * (t79 + (t83 ^ 2 + t86 ^ 2) * t75 ^ 2) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t73 ^ 2) / 0.2e1 + m(3) * (t77 ^ 2 + t79) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + qJ(2) * mrSges(3,3)) * t87 + (t61 * mrSges(5,1) - t62 * mrSges(5,2) + Ifges(5,3) * t78 / 0.2e1) * t78 + (t57 * mrSges(6,1) - t58 * mrSges(6,2) + Ifges(6,3) * t70 / 0.2e1) * t70 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t86 * mrSges(4,1) - t83 * mrSges(4,2)) * t75) * qJD(3) + (t77 * mrSges(3,2) + (mrSges(4,2) * t80 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t90 - t91) * t86) * t86 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t86) * qJD(1) + (Ifges(4,2) * t90 - t91) * t83) * t83) * qJD(1) + (t73 * mrSges(5,2) - t61 * mrSges(5,3) + Ifges(5,5) * t78 + Ifges(5,1) * t72 / 0.2e1) * t72 + (t59 * mrSges(6,2) - t57 * mrSges(6,3) + Ifges(6,5) * t70 + Ifges(6,1) * t65 / 0.2e1) * t65 - (-t73 * mrSges(5,1) + t62 * mrSges(5,3) + Ifges(5,4) * t72 + Ifges(5,6) * t78 - Ifges(5,2) * t71 / 0.2e1) * t71 + (-t59 * mrSges(6,1) + t58 * mrSges(6,3) + Ifges(6,4) * t65 + Ifges(6,6) * t70 + Ifges(6,2) * t64 / 0.2e1) * t64;
T = t1;
