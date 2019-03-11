% Return the minimum parameter vector for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRPRRR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t218 = (m(6) + m(7));
t209 = m(5) + t218;
t238 = (pkin(10) * t209);
t213 = sin(pkin(8));
t229 = mrSges(5,3) + t238;
t237 = t229 * t213;
t216 = cos(pkin(8));
t236 = t229 * t216;
t224 = (pkin(4) ^ 2);
t204 = (t224 * t218 + Ifges(5,2));
t226 = t204 + (2 * mrSges(5,3) + t238) * pkin(10);
t223 = (pkin(5) ^ 2);
t235 = (t223 * m(7) + Ifges(6,2));
t234 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t233 = pkin(3) ^ 2 * t209;
t231 = 2 * pkin(11) * mrSges(6,3) + t235;
t230 = pkin(12) * m(7) + mrSges(7,3);
t228 = pkin(11) * t218 + mrSges(6,3);
t221 = pkin(11) ^ 2;
t220 = pkin(12) ^ 2;
t217 = cos(pkin(7));
t215 = cos(pkin(14));
t214 = sin(pkin(7));
t212 = sin(pkin(14));
t1 = [m(2) + m(3); Ifges(3,3) + t217 ^ 2 * (t226 * t213 ^ 2 + Ifges(4,3) + t233) + (0.2e1 * (t212 * (-pkin(3) * t236 + Ifges(4,5)) + t215 * (t226 * t216 * t213 + Ifges(4,6))) * t217 + (t215 ^ 2 * (t226 * t216 ^ 2 + Ifges(4,2) + t233) + (0.2e1 * t215 * (pkin(3) * t237 + Ifges(4,4)) + (Ifges(4,1) + t226) * t212) * t212) * t214) * t214; mrSges(3,1); mrSges(3,2); pkin(3) * t209 + mrSges(4,1); mrSges(4,2) - t237; mrSges(4,3) + t236; m(4) + t209; t221 * t218 + Ifges(5,1) - t204 + t231; t228 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t221 + t224) * t218 + t231; pkin(4) * t218 + mrSges(5,1); mrSges(5,2) - t228; m(7) * t220 + Ifges(6,1) + t234 - t235; t230 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t220 + t223) * m(7) + t234; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t230; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
