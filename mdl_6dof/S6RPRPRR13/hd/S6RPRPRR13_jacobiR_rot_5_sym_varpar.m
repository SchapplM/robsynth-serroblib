% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR13_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:39
% EndTime: 2019-02-26 20:55:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (102->30), mult. (307->57), div. (0->0), fcn. (430->12), ass. (0->40)
t158 = cos(pkin(6));
t153 = sin(pkin(12));
t164 = cos(qJ(1));
t168 = t164 * t153;
t156 = cos(pkin(12));
t161 = sin(qJ(1));
t169 = t161 * t156;
t149 = t158 * t168 + t169;
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t167 = t164 * t156;
t170 = t161 * t153;
t148 = -t158 * t167 + t170;
t154 = sin(pkin(7));
t157 = cos(pkin(7));
t155 = sin(pkin(6));
t172 = t155 * t164;
t166 = t148 * t157 + t154 * t172;
t137 = t149 * t160 + t166 * t163;
t144 = -t148 * t154 + t157 * t172;
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t181 = -t137 * t159 + t144 * t162;
t180 = t137 * t162 + t144 * t159;
t177 = -t149 * t163 + t166 * t160;
t174 = t154 * t158;
t173 = t155 * t161;
t171 = t156 * t157;
t150 = -t158 * t169 - t168;
t165 = t150 * t157 + t154 * t173;
t151 = -t158 * t170 + t167;
t147 = -t155 * t156 * t154 + t158 * t157;
t146 = -t150 * t154 + t157 * t173;
t143 = t160 * t174 + (t153 * t163 + t160 * t171) * t155;
t142 = -t163 * t174 + (t153 * t160 - t163 * t171) * t155;
t141 = t151 * t163 + t165 * t160;
t140 = t151 * t160 - t165 * t163;
t136 = t140 * t159 + t146 * t162;
t135 = t140 * t162 - t146 * t159;
t1 = [t181, 0, t141 * t159, 0, t135, 0; t136, 0, -t177 * t159, 0, t180, 0; 0, 0, t143 * t159, 0, t142 * t162 - t147 * t159, 0; -t180, 0, t141 * t162, 0, -t136, 0; t135, 0, -t177 * t162, 0, t181, 0; 0, 0, t143 * t162, 0, -t142 * t159 - t147 * t162, 0; t177, 0, -t140, 0, 0, 0; t141, 0, -t137, 0, 0, 0; 0, 0, -t142, 0, 0, 0;];
JR_rot  = t1;
