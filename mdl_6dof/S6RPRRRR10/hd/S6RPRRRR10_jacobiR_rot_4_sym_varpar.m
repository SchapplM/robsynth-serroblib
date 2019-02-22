% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:04
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR10_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:04:18
% EndTime: 2019-02-22 11:04:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (105->30), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->40)
t160 = cos(pkin(6));
t155 = sin(pkin(13));
t166 = cos(qJ(1));
t170 = t166 * t155;
t158 = cos(pkin(13));
t163 = sin(qJ(1));
t171 = t163 * t158;
t150 = t160 * t170 + t171;
t162 = sin(qJ(3));
t165 = cos(qJ(3));
t169 = t166 * t158;
t172 = t163 * t155;
t149 = -t160 * t169 + t172;
t156 = sin(pkin(7));
t159 = cos(pkin(7));
t157 = sin(pkin(6));
t174 = t157 * t166;
t167 = t149 * t159 + t156 * t174;
t139 = -t150 * t165 + t167 * t162;
t144 = -t149 * t156 + t159 * t174;
t161 = sin(qJ(4));
t164 = cos(qJ(4));
t181 = t139 * t161 - t144 * t164;
t180 = t139 * t164 + t144 * t161;
t176 = t156 * t160;
t175 = t157 * t163;
t173 = t159 * t165;
t168 = t156 * t175;
t137 = -t150 * t162 - t167 * t165;
t152 = -t160 * t172 + t169;
t151 = -t160 * t171 - t170;
t148 = -t157 * t158 * t156 + t160 * t159;
t146 = -t151 * t156 + t159 * t175;
t143 = t162 * t176 + (t158 * t159 * t162 + t155 * t165) * t157;
t142 = t165 * t176 + (-t155 * t162 + t158 * t173) * t157;
t141 = t152 * t165 + (t151 * t159 + t168) * t162;
t140 = -t151 * t173 + t152 * t162 - t165 * t168;
t136 = t141 * t164 + t146 * t161;
t135 = -t141 * t161 + t146 * t164;
t1 = [t180, 0, -t140 * t164, t135, 0, 0; t136, 0, t137 * t164, t181, 0, 0; 0, 0, t142 * t164, -t143 * t161 + t148 * t164, 0, 0; -t181, 0, t140 * t161, -t136, 0, 0; t135, 0, -t137 * t161, t180, 0, 0; 0, 0, -t142 * t161, -t143 * t164 - t148 * t161, 0, 0; t137, 0, t141, 0, 0, 0; t140, 0, -t139, 0, 0, 0; 0, 0, t143, 0, 0, 0;];
JR_rot  = t1;
