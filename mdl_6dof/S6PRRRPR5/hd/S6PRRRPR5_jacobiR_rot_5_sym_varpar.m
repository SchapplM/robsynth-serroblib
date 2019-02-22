% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:56
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:56:06
% EndTime: 2019-02-22 09:56:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (130->39), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->46)
t163 = qJ(4) + pkin(13);
t161 = sin(t163);
t165 = sin(pkin(7));
t189 = t161 * t165;
t162 = cos(t163);
t188 = t162 * t165;
t166 = sin(pkin(6));
t187 = t165 * t166;
t169 = cos(pkin(6));
t186 = t165 * t169;
t168 = cos(pkin(7));
t185 = t166 * t168;
t170 = sin(qJ(3));
t184 = t168 * t170;
t172 = cos(qJ(3));
t183 = t168 * t172;
t171 = sin(qJ(2));
t182 = t169 * t171;
t173 = cos(qJ(2));
t181 = t169 * t173;
t180 = t170 * t171;
t179 = t170 * t173;
t178 = t171 * t172;
t177 = t172 * t173;
t176 = t171 * t187;
t164 = sin(pkin(12));
t167 = cos(pkin(12));
t156 = -t164 * t171 + t167 * t181;
t175 = t156 * t168 - t167 * t187;
t158 = -t164 * t181 - t167 * t171;
t174 = t158 * t168 + t164 * t187;
t159 = -t164 * t182 + t167 * t173;
t157 = t164 * t173 + t167 * t182;
t155 = t169 * t168 - t173 * t187;
t154 = (-t168 * t180 + t177) * t166;
t153 = -t158 * t165 + t164 * t185;
t152 = -t156 * t165 - t167 * t185;
t151 = t170 * t186 + (t168 * t179 + t178) * t166;
t150 = t172 * t186 + (t168 * t177 - t180) * t166;
t149 = t158 * t172 - t159 * t184;
t148 = t156 * t172 - t157 * t184;
t147 = t159 * t172 + t174 * t170;
t146 = -t159 * t170 + t174 * t172;
t145 = t157 * t172 + t175 * t170;
t144 = -t157 * t170 + t175 * t172;
t1 = [0, t149 * t162 + t159 * t189, t146 * t162, -t147 * t161 + t153 * t162, 0, 0; 0, t148 * t162 + t157 * t189, t144 * t162, -t145 * t161 + t152 * t162, 0, 0; 0, t154 * t162 + t161 * t176, t150 * t162, -t151 * t161 + t155 * t162, 0, 0; 0, -t149 * t161 + t159 * t188, -t146 * t161, -t147 * t162 - t153 * t161, 0, 0; 0, -t148 * t161 + t157 * t188, -t144 * t161, -t145 * t162 - t152 * t161, 0, 0; 0, -t154 * t161 + t162 * t176, -t150 * t161, -t151 * t162 - t155 * t161, 0, 0; 0, t158 * t170 + t159 * t183, t147, 0, 0, 0; 0, t156 * t170 + t157 * t183, t145, 0, 0, 0; 0 (t168 * t178 + t179) * t166, t151, 0, 0, 0;];
JR_rot  = t1;
