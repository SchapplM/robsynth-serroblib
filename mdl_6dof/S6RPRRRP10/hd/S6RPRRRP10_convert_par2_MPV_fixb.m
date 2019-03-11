% Return the minimum parameter vector for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRP10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t114 = (m(5) + m(6));
t112 = m(4) + t114;
t124 = (pkin(7) * t112);
t123 = (-Ifges(6,2) - Ifges(7,3));
t115 = (pkin(9) ^ 2);
t117 = (pkin(4) ^ 2);
t122 = (Ifges(5,2) + (t115 + t117) * m(6));
t121 = -pkin(9) * m(6) - mrSges(6,3);
t110 = (mrSges(5,3) - t121);
t120 = pkin(8) * t114 + t110;
t119 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(8) * t110 + t122 - t123;
t118 = (pkin(3) ^ 2);
t116 = pkin(8) ^ 2;
t1 = [t118 * t114 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + (2 * mrSges(4,3) + t124) * pkin(7); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t124; mrSges(3,3); m(3) + t112; Ifges(4,1) - Ifges(4,2) + (t116 - t118) * t114 + t119; t120 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t116 + t118) * t114 + t119; pkin(3) * t114 + mrSges(4,1); mrSges(4,2) - t120; t115 * m(6) + Ifges(5,1) - t122; Ifges(5,4); t121 * pkin(4) + Ifges(5,5); Ifges(5,6); t117 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t123; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
