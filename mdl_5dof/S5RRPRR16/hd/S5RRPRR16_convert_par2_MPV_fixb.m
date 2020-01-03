% Return the minimum parameter vector for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR16_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t132 = (m(3) * pkin(7));
t120 = (m(5) + m(6));
t131 = (pkin(8) * mrSges(5,3));
t124 = (pkin(4) ^ 2);
t130 = (t124 * m(6) + Ifges(5,2));
t129 = 2 * pkin(9) * mrSges(6,3) + Ifges(6,2);
t128 = pkin(9) * m(6) + mrSges(6,3);
t127 = -pkin(8) * t120 - mrSges(5,3);
t126 = (Ifges(3,2) + Ifges(4,3) + t130);
t125 = (pkin(3) ^ 2);
t123 = (pkin(8) ^ 2);
t122 = pkin(9) ^ 2;
t119 = sin(pkin(5));
t118 = 2 * t131;
t115 = (t123 + t125);
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (t115 * t120 + t118 + (2 * mrSges(3,3) + t132) * pkin(7) + t126) * t119 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t132) * t119; -2 * t131 + Ifges(3,1) + Ifges(4,2) + (-t115 + t125) * t120 - t126; Ifges(3,4) + Ifges(4,6); t127 * pkin(3) - Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,5); t123 * t120 + Ifges(4,1) + Ifges(3,3) + t118 + t130; mrSges(3,1); mrSges(3,2); pkin(3) * t120 + mrSges(4,1); mrSges(4,2) + t127; mrSges(4,3); m(4) + t120; m(6) * t122 + Ifges(5,1) + t129 - t130; t128 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t122 + t124) * m(6) + t129; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t128; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
