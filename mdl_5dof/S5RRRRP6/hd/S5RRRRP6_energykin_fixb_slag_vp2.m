% Calculate kinetic energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:45
% EndTime: 2019-12-31 21:52:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (356->82), mult. (752->120), div. (0->0), fcn. (486->6), ass. (0->28)
t93 = qJD(1) * (-pkin(7) - pkin(6));
t83 = sin(qJ(3));
t84 = sin(qJ(2));
t86 = cos(qJ(3));
t87 = cos(qJ(2));
t74 = (t83 * t84 - t86 * t87) * qJD(1);
t91 = pkin(6) * mrSges(3,3);
t75 = (t83 * t87 + t84 * t86) * qJD(1);
t80 = (-pkin(2) * t87 - pkin(1)) * qJD(1);
t64 = pkin(3) * t74 - pkin(8) * t75 + t80;
t78 = qJD(2) * pkin(2) + t84 * t93;
t79 = t87 * t93;
t69 = t83 * t78 - t86 * t79;
t81 = qJD(2) + qJD(3);
t67 = pkin(8) * t81 + t69;
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t60 = t82 * t64 + t85 * t67;
t59 = t85 * t64 - t67 * t82;
t68 = t78 * t86 + t83 * t79;
t66 = -pkin(3) * t81 - t68;
t73 = qJD(4) + t74;
t71 = t75 * t85 + t81 * t82;
t70 = -t75 * t82 + t81 * t85;
t61 = -pkin(4) * t70 + qJD(5) + t66;
t58 = qJ(5) * t70 + t60;
t57 = pkin(4) * t73 - qJ(5) * t71 + t59;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t80 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t61 ^ 2) / 0.2e1 + (t68 * mrSges(4,1) - t69 * mrSges(4,2) + Ifges(4,3) * t81 / 0.2e1) * t81 + (t80 * mrSges(4,2) - t68 * mrSges(4,3) + Ifges(4,5) * t81 + Ifges(4,1) * t75 / 0.2e1) * t75 - (-t80 * mrSges(4,1) + t69 * mrSges(4,3) + Ifges(4,4) * t75 + Ifges(4,6) * t81 - Ifges(4,2) * t74 / 0.2e1) * t74 + (t59 * mrSges(5,1) + t57 * mrSges(6,1) - t60 * mrSges(5,2) - t58 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t73) * t73 + (t66 * mrSges(5,2) + t61 * mrSges(6,2) - t59 * mrSges(5,3) - t57 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t71 + (Ifges(5,5) + Ifges(6,5)) * t73) * t71 + (-t66 * mrSges(5,1) - t61 * mrSges(6,1) + t60 * mrSges(5,3) + t58 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t70 + (Ifges(5,6) + Ifges(6,6)) * t73 + (Ifges(5,4) + Ifges(6,4)) * t71) * t70 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t84 ^ 2 + t87 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t91 + Ifges(3,2) / 0.2e1) * t87) * t87 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t87 + (t91 + Ifges(3,1) / 0.2e1) * t84) * t84) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t87 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t84) * qJD(2)) * qJD(1);
T = t1;
