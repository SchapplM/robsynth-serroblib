% Return the minimum parameter vector for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRRRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t147 = (m(5) + m(6));
t143 = m(4) + t147;
t158 = (t143 * pkin(8));
t151 = (pkin(4) ^ 2);
t157 = (t151 * m(6) + Ifges(5,2));
t156 = 2 * pkin(10) * mrSges(6,3) + Ifges(6,2);
t155 = 2 * pkin(9) * mrSges(5,3) + t157;
t154 = pkin(10) * m(6) + mrSges(6,3);
t153 = pkin(9) * t147 + mrSges(5,3);
t152 = (pkin(3) ^ 2);
t150 = pkin(9) ^ 2;
t149 = pkin(10) ^ 2;
t146 = sin(pkin(6));
t1 = [m(2) + m(3) + t143; pkin(2) ^ 2 * t143 + Ifges(3,3) + (t152 * t147 + Ifges(4,2) + (2 * mrSges(4,3) + t158) * pkin(8)) * t146 ^ 2; pkin(2) * t143 + mrSges(3,1); mrSges(3,2) + (-mrSges(4,3) - t158) * t146; Ifges(4,1) - Ifges(4,2) + (t150 - t152) * t147 + t155; t153 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t150 + t152) * t147 + t155; pkin(3) * t147 + mrSges(4,1); mrSges(4,2) - t153; m(6) * t149 + Ifges(5,1) + t156 - t157; t154 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t149 + t151) * m(6) + t156; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t154; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
