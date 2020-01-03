% Calculate kinetic energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:18
% EndTime: 2020-01-03 11:57:19
% DurationCPUTime: 0.19s
% Computational Cost: add. (208->51), mult. (327->86), div. (0->0), fcn. (160->8), ass. (0->26)
t78 = qJD(1) + qJD(2);
t86 = cos(qJ(2));
t91 = qJD(1) * pkin(1);
t72 = t78 * pkin(2) + t86 * t91;
t80 = sin(pkin(8));
t82 = cos(pkin(8));
t84 = sin(qJ(2));
t90 = t84 * t91;
t70 = t80 * t72 + t82 * t90;
t68 = qJ(4) * t78 + t70;
t79 = sin(pkin(9));
t81 = cos(pkin(9));
t64 = -t81 * qJD(3) + t68 * t79;
t93 = t64 ^ 2;
t92 = t78 * t81;
t69 = t72 * t82 - t80 * t90;
t89 = qJD(4) - t69;
t85 = cos(qJ(5));
t83 = sin(qJ(5));
t73 = qJD(5) - t92;
t67 = -pkin(3) * t78 + t89;
t66 = qJD(3) * t79 + t68 * t81;
t63 = (-pkin(4) * t81 - pkin(7) * t79 - pkin(3)) * t78 + t89;
t62 = t83 * t63 + t85 * t66;
t61 = t85 * t63 - t83 * t66;
t1 = m(5) * (t66 ^ 2 + t67 ^ 2 + t93) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t93) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t84 ^ 2 + t86 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t73 / 0.2e1) * t73 + (t69 * mrSges(4,1) - t70 * mrSges(4,2) + (mrSges(3,1) * t86 - mrSges(3,2) * t84) * t91 + (-t67 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,2) * t92 / 0.2e1) * t81 + (t67 * mrSges(5,2) + (mrSges(6,1) * t83 + mrSges(6,2) * t85 + mrSges(5,3)) * t64 + (-t61 * t85 - t62 * t83) * mrSges(6,3) + t73 * (Ifges(6,5) * t85 - Ifges(6,6) * t83)) * t79 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t81 + (Ifges(6,1) * t85 ^ 2 / 0.2e1 + Ifges(5,1) / 0.2e1 + (-Ifges(6,4) * t85 + Ifges(6,2) * t83 / 0.2e1) * t83) * t79) * t79) * t78) * t78;
T = t1;
