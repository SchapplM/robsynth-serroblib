% Calculate kinetic energy for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:37
% EndTime: 2019-12-31 21:50:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (292->67), mult. (409->103), div. (0->0), fcn. (210->6), ass. (0->27)
t74 = qJD(1) + qJD(2);
t88 = t74 / 0.2e1;
t87 = cos(qJ(4));
t77 = sin(qJ(2));
t84 = pkin(1) * qJD(1);
t70 = pkin(7) * t74 + t77 * t84;
t86 = t70 * mrSges(4,3);
t78 = cos(qJ(3));
t85 = t74 * t78;
t76 = sin(qJ(3));
t82 = pkin(8) * t74 + t70;
t64 = qJD(3) * pkin(3) - t82 * t76;
t65 = t82 * t78;
t75 = sin(qJ(4));
t61 = t75 * t64 + t87 * t65;
t79 = cos(qJ(2));
t83 = t79 * t84;
t60 = t87 * t64 - t75 * t65;
t68 = -t83 + (-pkin(3) * t78 - pkin(2)) * t74;
t73 = qJD(3) + qJD(4);
t71 = -pkin(2) * t74 - t83;
t67 = (t75 * t78 + t87 * t76) * t74;
t66 = t74 * t75 * t76 - t87 * t85;
t59 = pkin(4) * t66 - qJ(5) * t67 + t68;
t58 = qJ(5) * t73 + t61;
t57 = -t73 * pkin(4) + qJD(5) - t60;
t1 = m(4) * (t71 ^ 2 + (t76 ^ 2 + t78 ^ 2) * t70 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t77 ^ 2 + t79 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-t76 * mrSges(4,1) - t78 * mrSges(4,2)) * t70) * qJD(3) + (Ifges(3,3) * t88 + (mrSges(3,1) * t79 - mrSges(3,2) * t77) * t84 + (-t71 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t88 + t86) * t78) * t78 + (Ifges(4,4) * t85 + t71 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t88 + t86) * t76) * t76) * t74 + (t60 * mrSges(5,1) - t57 * mrSges(6,1) - t61 * mrSges(5,2) + t58 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t73) * t73 + (t68 * mrSges(5,2) + t57 * mrSges(6,2) - t60 * mrSges(5,3) - t59 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t67 + (Ifges(6,4) + Ifges(5,5)) * t73) * t67 + (t68 * mrSges(5,1) + t59 * mrSges(6,1) - t58 * mrSges(6,2) - t61 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t66 + (-Ifges(5,6) + Ifges(6,6)) * t73 + (-Ifges(5,4) + Ifges(6,5)) * t67) * t66;
T = t1;
