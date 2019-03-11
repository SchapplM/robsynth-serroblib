% Return the minimum parameter vector for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRP14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t163 = (m(3) * pkin(8));
t150 = (m(5) + m(6));
t162 = (pkin(9) * mrSges(5,3));
t161 = (-Ifges(6,2) - Ifges(7,3));
t154 = (pkin(4) ^ 2);
t160 = (t154 * m(6) + Ifges(5,2));
t159 = 2 * pkin(10) * mrSges(6,3) - t161;
t158 = pkin(10) * m(6) + mrSges(6,3);
t157 = -pkin(9) * t150 - mrSges(5,3);
t156 = (Ifges(3,2) + Ifges(4,3) + t160);
t155 = (pkin(3) ^ 2);
t153 = (pkin(9) ^ 2);
t152 = pkin(10) ^ 2;
t149 = sin(pkin(6));
t148 = 2 * t162;
t145 = (t153 + t155);
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (t145 * t150 + t148 + (2 * mrSges(3,3) + t163) * pkin(8) + t156) * t149 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t163) * t149; -2 * t162 + Ifges(3,1) + Ifges(4,2) + (-t145 + t155) * t150 - t156; Ifges(3,4) + Ifges(4,6); t157 * pkin(3) - Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,5); t153 * t150 + Ifges(4,1) + Ifges(3,3) + t148 + t160; mrSges(3,1); mrSges(3,2); pkin(3) * t150 + mrSges(4,1); mrSges(4,2) + t157; mrSges(4,3); m(4) + t150; t152 * m(6) + Ifges(5,1) + t159 - t160; t158 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t152 + t154) * m(6) + t159; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t158; Ifges(6,1) + Ifges(7,1) + t161; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
