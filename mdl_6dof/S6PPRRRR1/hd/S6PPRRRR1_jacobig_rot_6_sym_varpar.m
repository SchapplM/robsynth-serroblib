% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (49->21), mult. (129->46), div. (0->0), fcn. (184->12), ass. (0->32)
t166 = sin(pkin(12));
t172 = cos(pkin(6));
t182 = t166 * t172;
t167 = sin(pkin(7));
t168 = sin(pkin(6));
t181 = t167 * t168;
t180 = t167 * t172;
t171 = cos(pkin(7));
t179 = t168 * t171;
t169 = cos(pkin(13));
t178 = t169 * t171;
t170 = cos(pkin(12));
t177 = t170 * t172;
t165 = sin(pkin(13));
t158 = -t165 * t166 + t169 * t177;
t176 = -t158 * t171 + t170 * t181;
t160 = -t165 * t170 - t169 * t182;
t175 = t160 * t171 + t166 * t181;
t174 = cos(qJ(3));
t173 = sin(qJ(3));
t164 = qJ(4) + qJ(5);
t163 = cos(t164);
t162 = sin(t164);
t161 = -t165 * t182 + t169 * t170;
t159 = t165 * t177 + t166 * t169;
t157 = -t169 * t181 + t171 * t172;
t156 = -t160 * t167 + t166 * t179;
t155 = -t158 * t167 - t170 * t179;
t154 = -t174 * t180 + (t165 * t173 - t174 * t178) * t168;
t153 = t161 * t173 - t174 * t175;
t152 = t159 * t173 + t174 * t176;
t1 = [0, 0, t156, t153, t153 (t161 * t174 + t173 * t175) * t162 - t156 * t163; 0, 0, t155, t152, t152 (t159 * t174 - t173 * t176) * t162 - t155 * t163; 0, 0, t157, t154, t154 (t173 * t180 + (t165 * t174 + t173 * t178) * t168) * t162 - t157 * t163;];
Jg_rot  = t1;
