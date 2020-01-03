% Calculate kinetic energy for
% S5RRRRP3
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:04
% EndTime: 2019-12-31 21:49:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (205->53), mult. (265->83), div. (0->0), fcn. (104->6), ass. (0->20)
t69 = qJD(1) + qJD(2);
t75 = cos(qJ(2));
t79 = pkin(1) * qJD(1);
t65 = pkin(2) * t69 + t75 * t79;
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t72 = sin(qJ(2));
t78 = t72 * t79;
t63 = t71 * t65 + t74 * t78;
t68 = qJD(3) + t69;
t61 = pkin(8) * t68 + t63;
t80 = t61 * mrSges(5,3);
t62 = t65 * t74 - t71 * t78;
t73 = cos(qJ(4));
t70 = sin(qJ(4));
t60 = -pkin(3) * t68 - t62;
t58 = qJD(4) * qJ(5) + t61 * t73;
t57 = -qJD(4) * pkin(4) + t61 * t70 + qJD(5);
t56 = (-pkin(4) * t73 - qJ(5) * t70 - pkin(3)) * t68 - t62;
t1 = m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + (t70 ^ 2 + t73 ^ 2) * t61 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * t69 / 0.2e1 + (mrSges(3,1) * t75 - mrSges(3,2) * t72) * t79) * t69 + (-t57 * mrSges(6,1) + t58 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4) + (-mrSges(5,1) * t70 - mrSges(5,2) * t73) * t61) * qJD(4) + (-t63 * mrSges(4,2) + t62 * mrSges(4,1) + Ifges(4,3) * t68 / 0.2e1 + (-t60 * mrSges(5,1) - t56 * mrSges(6,1) + t58 * mrSges(6,2) + (t80 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t68) * t73 + (Ifges(5,6) - Ifges(6,6)) * qJD(4)) * t73 + (t60 * mrSges(5,2) - t56 * mrSges(6,3) + t57 * mrSges(6,2) + (t80 + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t68) * t70 + (Ifges(5,4) - Ifges(6,5)) * t68 * t73 + (Ifges(6,4) + Ifges(5,5)) * qJD(4)) * t70) * t68;
T = t1;
