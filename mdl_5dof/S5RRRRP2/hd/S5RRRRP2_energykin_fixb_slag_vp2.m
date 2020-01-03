% Calculate kinetic energy for
% S5RRRRP2
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:01
% EndTime: 2020-01-03 12:11:01
% DurationCPUTime: 0.21s
% Computational Cost: add. (296->67), mult. (417->103), div. (0->0), fcn. (218->6), ass. (0->26)
t74 = qJD(1) + qJD(2);
t87 = t74 / 0.2e1;
t77 = sin(qJ(2));
t85 = qJD(1) * pkin(1);
t71 = pkin(7) * t74 + t77 * t85;
t86 = t71 * mrSges(4,3);
t76 = sin(qJ(3));
t83 = pkin(8) * t74 + t71;
t65 = qJD(3) * pkin(3) - t83 * t76;
t79 = cos(qJ(3));
t66 = t83 * t79;
t75 = sin(qJ(4));
t78 = cos(qJ(4));
t60 = t75 * t65 + t78 * t66;
t80 = cos(qJ(2));
t84 = t80 * t85;
t59 = t78 * t65 - t66 * t75;
t69 = -t84 + (-pkin(3) * t79 - pkin(2)) * t74;
t73 = qJD(3) + qJD(4);
t72 = -pkin(2) * t74 - t84;
t68 = (t75 * t79 + t76 * t78) * t74;
t67 = (-t75 * t76 + t78 * t79) * t74;
t61 = -pkin(4) * t67 + qJD(5) + t69;
t58 = qJ(5) * t67 + t60;
t57 = pkin(4) * t73 - qJ(5) * t68 + t59;
t1 = m(4) * (t72 ^ 2 + (t76 ^ 2 + t79 ^ 2) * t71 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t61 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t77 ^ 2 + t80 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-t76 * mrSges(4,1) - t79 * mrSges(4,2)) * t71) * qJD(3) + (Ifges(3,3) * t87 + (mrSges(3,1) * t80 - mrSges(3,2) * t77) * t85 + (-t72 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t87 + t86) * t79) * t79 + (Ifges(4,4) * t74 * t79 + t72 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t87 + t86) * t76) * t76) * t74 + (t59 * mrSges(5,1) + t57 * mrSges(6,1) - t60 * mrSges(5,2) - t58 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t73) * t73 + (t69 * mrSges(5,2) + t61 * mrSges(6,2) - t59 * mrSges(5,3) - t57 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t68 + (Ifges(5,5) + Ifges(6,5)) * t73) * t68 + (-t69 * mrSges(5,1) - t61 * mrSges(6,1) + t60 * mrSges(5,3) + t58 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t67 + (Ifges(5,6) + Ifges(6,6)) * t73 + (Ifges(5,4) + Ifges(6,4)) * t68) * t67;
T = t1;
