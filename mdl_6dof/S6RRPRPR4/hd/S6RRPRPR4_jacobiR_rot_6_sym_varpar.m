% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:24:31
% EndTime: 2019-02-22 11:24:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (205->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->38)
t177 = qJ(4) + pkin(12);
t175 = sin(t177);
t176 = cos(t177);
t181 = cos(pkin(6));
t178 = sin(pkin(11));
t180 = cos(pkin(11));
t183 = sin(qJ(2));
t186 = cos(qJ(2));
t189 = t186 * t178 + t183 * t180;
t169 = t189 * t181;
t170 = t183 * t178 - t186 * t180;
t184 = sin(qJ(1));
t187 = cos(qJ(1));
t191 = t187 * t169 - t184 * t170;
t179 = sin(pkin(6));
t194 = t179 * t187;
t155 = t175 * t194 - t176 * t191;
t188 = t170 * t181;
t159 = -t184 * t189 - t187 * t188;
t182 = sin(qJ(6));
t185 = cos(qJ(6));
t201 = t155 * t182 - t159 * t185;
t200 = t155 * t185 + t159 * t182;
t197 = t176 * t182;
t196 = t176 * t185;
t195 = t179 * t184;
t190 = -t184 * t169 - t187 * t170;
t153 = -t175 * t191 - t176 * t194;
t168 = t189 * t179;
t167 = t170 * t179;
t165 = t168 * t176 + t181 * t175;
t164 = -t168 * t175 + t181 * t176;
t162 = t184 * t188 - t187 * t189;
t157 = t175 * t195 + t176 * t190;
t156 = t175 * t190 - t176 * t195;
t152 = t157 * t185 - t162 * t182;
t151 = -t157 * t182 - t162 * t185;
t1 = [t200, t162 * t196 + t182 * t190, 0, -t156 * t185, 0, t151; t152, t159 * t196 + t182 * t191, 0, t153 * t185, 0, t201; 0, -t167 * t196 + t168 * t182, 0, t164 * t185, 0, -t165 * t182 + t167 * t185; -t201, -t162 * t197 + t185 * t190, 0, t156 * t182, 0, -t152; t151, -t159 * t197 + t185 * t191, 0, -t153 * t182, 0, t200; 0, t167 * t197 + t168 * t185, 0, -t164 * t182, 0, -t165 * t185 - t167 * t182; t153, t162 * t175, 0, t157, 0, 0; t156, t159 * t175, 0, -t155, 0, 0; 0, -t167 * t175, 0, t165, 0, 0;];
JR_rot  = t1;
