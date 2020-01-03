% Calculate kinetic energy for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:40
% EndTime: 2019-12-31 18:17:40
% DurationCPUTime: 0.15s
% Computational Cost: add. (120->44), mult. (212->68), div. (0->0), fcn. (80->6), ass. (0->21)
t69 = cos(pkin(8));
t64 = (pkin(1) * t69 + pkin(2)) * qJD(1);
t71 = sin(qJ(3));
t73 = cos(qJ(3));
t68 = sin(pkin(8));
t78 = pkin(1) * qJD(1) * t68;
t63 = t64 * t71 + t73 * t78;
t67 = qJD(1) + qJD(3);
t60 = qJ(4) * t67 + t63;
t80 = t60 ^ 2;
t79 = m(3) / 0.2e1;
t62 = t73 * t64 - t71 * t78;
t77 = qJD(4) - t62;
t74 = qJD(2) ^ 2;
t72 = cos(qJ(5));
t70 = sin(qJ(5));
t59 = -pkin(3) * t67 + t77;
t58 = (-pkin(3) - pkin(7)) * t67 + t77;
t57 = qJD(2) * t72 + t58 * t70;
t56 = -qJD(2) * t70 + t58 * t72;
t1 = t74 * t79 + m(4) * (t62 ^ 2 + t63 ^ 2 + t74) / 0.2e1 + m(5) * (t59 ^ 2 + t74 + t80) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t80) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t69 * mrSges(3,1) - t68 * mrSges(3,2) + (t68 ^ 2 + t69 ^ 2) * t79 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (t62 * mrSges(4,1) - t63 * mrSges(4,2) + t59 * mrSges(5,2) + t60 * mrSges(5,3) + (t60 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t72 + (t60 * mrSges(6,1) - t57 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t70 + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) * t72 ^ 2 / 0.2e1 + (-Ifges(6,4) * t72 + Ifges(6,2) * t70 / 0.2e1) * t70) * t67) * t67;
T = t1;
