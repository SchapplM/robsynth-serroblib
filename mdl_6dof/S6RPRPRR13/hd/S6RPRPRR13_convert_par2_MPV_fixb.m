% Return the minimum parameter vector for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR13_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t220 = (m(4) * pkin(9));
t195 = sin(pkin(7));
t212 = mrSges(4,3) + t220;
t219 = t195 * t212;
t198 = cos(pkin(7));
t218 = t198 * t212;
t200 = (m(6) + m(7));
t217 = (pkin(10) * mrSges(6,3));
t216 = pkin(2) ^ 2 * m(4);
t205 = (pkin(5) ^ 2);
t215 = (t205 * m(7) + Ifges(6,2));
t214 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t211 = pkin(11) * m(7) + mrSges(7,3);
t210 = -pkin(10) * t200 - mrSges(6,3);
t203 = (pkin(10) ^ 2);
t206 = (pkin(4) ^ 2);
t209 = (Ifges(4,2) + Ifges(5,3) + (t203 + t206) * t200 + t215);
t193 = 2 * t217;
t208 = t193 + t209 + (2 * mrSges(4,3) + t220) * pkin(9);
t202 = pkin(11) ^ 2;
t199 = cos(pkin(6));
t197 = cos(pkin(12));
t196 = sin(pkin(6));
t194 = sin(pkin(12));
t1 = [Ifges(2,3) + t199 ^ 2 * (t195 ^ 2 * t208 + Ifges(3,3) + t216) + (0.2e1 * (t194 * (-pkin(2) * t218 + Ifges(3,5)) + t197 * (t195 * t198 * t208 + Ifges(3,6))) * t199 + (t197 ^ 2 * (t198 ^ 2 * t208 + Ifges(3,2) + t216) + (0.2e1 * t197 * (pkin(2) * t219 + Ifges(3,4)) + (Ifges(3,1) + t208) * t194) * t194) * t196) * t196; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t219; mrSges(3,3) + t218; m(3) + m(4); t206 * t200 + Ifges(4,1) + Ifges(5,2) - t209 - 2 * t217; Ifges(4,4) + Ifges(5,6); pkin(4) * t210 - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5); t203 * t200 + Ifges(5,1) + Ifges(4,3) + t193 + t215; mrSges(4,1); mrSges(4,2); pkin(4) * t200 + mrSges(5,1); mrSges(5,2) + t210; mrSges(5,3); m(5) + t200; m(7) * t202 + Ifges(6,1) + t214 - t215; pkin(5) * t211 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t202 + t205) * m(7) + t214; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t211; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
