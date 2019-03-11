% Return the minimum parameter vector for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRP11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t134 = -pkin(9) * m(6) - mrSges(6,3);
t122 = (mrSges(5,3) - t134);
t127 = (pkin(9) ^ 2);
t129 = (pkin(4) ^ 2);
t135 = (Ifges(5,2) + (t127 + t129) * m(6));
t136 = (-Ifges(6,2) - Ifges(7,2));
t131 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(8) * t122 + t135 - t136;
t140 = -t131 - Ifges(3,2) - Ifges(4,3);
t126 = (m(5) + m(6));
t133 = -pkin(8) * t126 - t122;
t130 = (pkin(3) ^ 2);
t128 = pkin(8) ^ 2;
t124 = t128 + t130;
t1 = [Ifges(2,3) + t124 * t126 + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * m(3) - t140; m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) + Ifges(4,2) + (-t124 + t130) * t126 + t140; Ifges(3,4) + Ifges(4,6); t133 * pkin(3) - Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,5); t128 * t126 + Ifges(4,1) + Ifges(3,3) + t131; mrSges(3,1); mrSges(3,2); pkin(3) * t126 + mrSges(4,1); mrSges(4,2) + t133; mrSges(4,3); m(4) + t126; t127 * m(6) + Ifges(5,1) - t135; Ifges(5,4); t134 * pkin(4) + Ifges(5,5); Ifges(5,6); t129 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t136; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
