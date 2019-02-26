% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:34
% EndTime: 2019-02-26 22:29:35
% DurationCPUTime: 0.15s
% Computational Cost: add. (74->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t150 = cos(pkin(6));
t153 = sin(qJ(2));
t158 = cos(qJ(1));
t161 = t158 * t153;
t154 = sin(qJ(1));
t157 = cos(qJ(2));
t164 = t154 * t157;
t144 = t150 * t161 + t164;
t152 = sin(qJ(3));
t156 = cos(qJ(3));
t149 = sin(pkin(6));
t167 = t149 * t158;
t137 = t144 * t156 - t152 * t167;
t160 = t158 * t157;
t165 = t154 * t153;
t143 = -t150 * t160 + t165;
t151 = sin(qJ(4));
t155 = cos(qJ(4));
t174 = -t137 * t151 + t143 * t155;
t173 = t137 * t155 + t143 * t151;
t170 = t149 * t152;
t169 = t149 * t156;
t168 = t149 * t157;
t166 = t151 * t156;
t163 = t155 * t156;
t162 = t156 * t157;
t159 = t144 * t152 + t156 * t167;
t146 = -t150 * t165 + t160;
t145 = t150 * t164 + t161;
t142 = t150 * t152 + t153 * t169;
t141 = t150 * t156 - t153 * t170;
t140 = t146 * t156 + t154 * t170;
t139 = -t146 * t152 + t154 * t169;
t135 = t140 * t155 + t145 * t151;
t134 = t140 * t151 - t145 * t155;
t1 = [-t173, -t145 * t163 + t146 * t151, t139 * t155, -t134, 0, 0; t135, -t143 * t163 + t144 * t151, -t159 * t155, t174, 0, 0; 0 (t151 * t153 + t155 * t162) * t149, t141 * t155, -t142 * t151 - t155 * t168, 0, 0; t174, -t145 * t166 - t146 * t155, t139 * t151, t135, 0, 0; t134, -t143 * t166 - t144 * t155, -t159 * t151, t173, 0, 0; 0 (t151 * t162 - t153 * t155) * t149, t141 * t151, t142 * t155 - t151 * t168, 0, 0; t159, t145 * t152, -t140, 0, 0, 0; t139, t143 * t152, -t137, 0, 0, 0; 0, -t152 * t168, -t142, 0, 0, 0;];
JR_rot  = t1;
